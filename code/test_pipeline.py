#!/usr/bin/env python3
"""
code/test_pipeline.py

Synthetic and adversarial unit tests for canonical_classifier.py and
run_saddle_replication.py.
"""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import canonical_classifier as cc
from run_saddle_replication import extract_strict_amplicon, require_bowtie2, write_pair_table_tsv


class TestCanonicalClassifier(unittest.TestCase):
    def test_adapter_trimming_exact_and_missing(self):
        adp1 = cc.ORIG_ADAPTER_R1
        adp2 = cc.ORIG_ADAPTER_R2
        s1 = "ACGTACGTACGTACGT" + adp1 + "TTTT"
        q1 = "I" * len(s1)
        s2 = "TGCATGCATGCATGCA" + adp2 + "AAAA"
        q2 = "I" * len(s2)

        t1, tq1, t2, tq2, st = cc.trim_adapter_pair(s1, q1, s2, q2, adp1, adp2)
        self.assertEqual(st, "ok")
        self.assertEqual(t1, "ACGTACGTACGTACGT")
        self.assertEqual(t2, "TGCATGCATGCATGCA")
        self.assertEqual(len(tq1), len(t1))
        self.assertEqual(len(tq2), len(t2))

        _, _, _, _, st_both = cc.trim_adapter_pair("AAAA", "IIII", "TTTT", "IIII", adp1, adp2)
        self.assertEqual(st_both, "no_adapter_both")
        _, _, _, _, st_r1 = cc.trim_adapter_pair("AAAA", "IIII", s2, q2, adp1, adp2)
        self.assertEqual(st_r1, "no_adapter_r1_only")
        _, _, _, _, st_r2 = cc.trim_adapter_pair(s1, q1, "TTTT", "IIII", adp1, adp2)
        self.assertEqual(st_r2, "no_adapter_r2_only")

    def test_qc_filters(self):
        self.assertEqual(cc.qc_check_pair("ACGTACGTAC", "I" * 10, "ACGTACGTACG", "I" * 11), "qc_short_len")
        self.assertEqual(cc.qc_check_pair("ACGTNCGTACG", "I" * 11, "ACGTACGTACG", "I" * 11), "qc_ambiguous_n")
        self.assertEqual(cc.qc_check_pair("ACGTACGTACG", "IIIIIII!!!!", "ACGTACGTACG", "I" * 11), "qc_low_quality")
        self.assertEqual(cc.qc_check_pair("ACGTACGTACG", "IIIIIIII!!!", "ACGTACGTACG", "I" * 11), "ok")

    def test_sw_nuc44_scoring_substitutions_gaps_and_boundaries(self):
        p = "gatctgctgtttcactttgggtt"  # 23 bp, perfect score = 115
        self.assertEqual(cc.sw_score_nuc44(p, p), 115)
        p_sub1 = p[:11] + "a" + p[12:]
        self.assertEqual(cc.sw_score_nuc44(p, p_sub1), 106)
        p_del1 = p[:11] + p[12:]
        self.assertEqual(cc.sw_score_nuc44(p, p_del1), 102)

        primers = [p]
        self.assertEqual(cc.match_primer_sw_nuc44(p + "acgtacgt", [0], primers), 0)
        self.assertIsNone(cc.match_primer_sw_nuc44(p[:14], [0], primers))
        self.assertEqual(cc.match_primer_prefix(p[:14], [0], primers), 0)

        p20 = "acgtacgtacgtacgtacgt"
        q_border = p20[:12] + "tttttttt"
        self.assertEqual(cc.sw_score_nuc44(p20, q_border), 60)
        self.assertIsNone(cc.match_primer_sw_nuc44(q_border, [0], [p20], threshold=0.6))
        q_pass = p20[:13] + "ttttttt"
        self.assertEqual(cc.match_primer_sw_nuc44(q_pass, [0], [p20], threshold=0.6), 0)

    def test_dimer_threshold(self):
        self.assertTrue(cc.is_dimer_insert(37, 20, 20))
        self.assertFalse(cc.is_dimer_insert(38, 20, 20))

    def test_edlib_recovers_frameshift_indel(self):
        ref = "gatctgctgtttcactttgggttacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt"
        read_del = ref[:25] + ref[26:]
        self.assertLess(cc.ungapped_positional_identity(read_del, ref), 0.70)
        self.assertGreater(cc.edlib_shw_identity(read_del, ref), 0.95)

    def test_adversarial_fastq_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            r1 = d / "ok_1.fastq"
            r2 = d / "ok_2.fastq"
            cc.write_utf8_lf(r1, "@readA/1 desc\nACGTACGT\n+\nIIIIIIII\n")
            cc.write_utf8_lf(r2, "@readA/2 desc\nTGCATGCA\n+\nIIIIIIII\n")
            recs = list(cc.iter_validated_fastq_pairs(r1, r2))
            self.assertEqual(len(recs), 1)
            self.assertEqual(recs[0][0], "@readA")

            cc.write_utf8_lf(r2, "@readB/2\nTGCATGCA\n+\nIIIIIIII\n")
            with self.assertRaises(ValueError):
                list(cc.iter_validated_fastq_pairs(r1, r2))

            cc.write_utf8_lf(r2, "@readA/2\nTGCATGCA\n+\n")
            with self.assertRaises(ValueError):
                list(cc.iter_validated_fastq_pairs(r1, r2))

            cc.write_utf8_lf(r2, "@readA/2\nTGCATGCA\n+\nIIII\n")
            with self.assertRaises(ValueError):
                list(cc.iter_validated_fastq_pairs(r1, r2))

            cc.write_utf8_lf(r1, "@r1/1\nACGT\n+\nIIII\n@r2/1\nACGT\n+\nIIII\n")
            cc.write_utf8_lf(r2, "@r1/2\nACGT\n+\nIIII\n")
            with self.assertRaises(ValueError):
                list(cc.iter_validated_fastq_pairs(r1, r2))

    def test_strict_amplicon_extraction_no_silent_fallback(self):
        fp = "gatctgctgtttcactttgggtt"
        rp = "cgtcaggctgtatttcttccac"
        rp_rc = cc.revcomp(rp)
        ctx = "aaaa" + fp + "cccccc" + rp_rc + "tttt"
        amp = extract_strict_amplicon("fp1", fp, rp, ctx)
        self.assertEqual(amp, fp + "cccccc" + rp_rc)
        with self.assertRaises(ValueError):
            extract_strict_amplicon("fp1", fp, rp, "aaaa" + fp + "cccccc")

    def test_directional_matrix_export_and_fail_fast_bowtie2(self):
        mat = np.array([[3, 7], [2, 5]], dtype=np.int64)
        primers = [("p0", "acgt"), ("p1", "tgca")]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "m.tsv.gz"
            write_pair_table_tsv(mat, primers, out, include_directional=True, compress_gz=True)
            with gzip.open(out, "rt", encoding="utf-8") as f:
                rows = [line.strip().split("\t") for line in f]
            self.assertEqual(rows[0], ["pair_index", "primer_1_id", "primer_2_id", "reads_p1_p2", "reads_p2_p1", "dimer_reads"])
            # 3 upper-triangular entries: (0,0)->3+0=3, (0,1)->7+2=9, (1,1)->5+0=5
            self.assertEqual(rows[1][3:], ["3", "0", "3"])
            self.assertEqual(rows[2][3:], ["7", "2", "9"])
            self.assertEqual(rows[3][3:], ["5", "0", "5"])

        with mock.patch("run_saddle_replication.get_bowtie2_paths", return_value=(None, None, None)):
            with self.assertRaises(RuntimeError):
                require_bowtie2(skip_bowtie2=False)
            run_bt2, _, _, _ = require_bowtie2(skip_bowtie2=True)
            self.assertFalse(run_bt2)

    def test_cross_target_and_tied_match_edge_cases(self):
        amp0 = "gatctgctgtttcactttgggttacgtacgtacgtacgtacgt"
        amp1 = "ttggccaaagggccctttaaacccgggaaatttcccgggaaat"
        amp0_rc = cc.revcomp(amp0)
        amp1_rc = cc.revcomp(amp1)

        # (a) Reverse-strand exact match to the SAME assigned amplicon (amp0) must NOT count as non-cognate
        r_same_rc = amp0_rc[:30]
        res_rc = cc.classify_best_target_edlib(r_same_rc, [amp0, amp1], [amp0_rc, amp1_rc], assigned_amp_idx=0)
        self.assertEqual(res_rc["best_amp_1based"], 1)
        self.assertEqual(res_rc["best_strand"], "-")
        self.assertFalse(res_rc["strictly_better_noncognate"])
        self.assertFalse(res_rc["tied_best_noncognate"])

        # (b) Strictly better match to a DIFFERENT amplicon (amp1) when assigned to amp0
        r_other = amp1[:30]
        res_other = cc.classify_best_target_edlib(r_other, [amp0, amp1], [amp0_rc, amp1_rc], assigned_amp_idx=0)
        self.assertEqual(res_other["best_amp_1based"], 2)
        self.assertTrue(res_other["strictly_better_noncognate"])
        self.assertFalse(res_other["tied_best_noncognate"])

        # (c) Exact tie between assigned amp0 and non-cognate amp1_twin
        res_tie = cc.classify_best_target_edlib(amp0[:30], [amp0, amp0], [amp0_rc, amp0_rc], assigned_amp_idx=0)
        self.assertFalse(res_tie["strictly_better_noncognate"])
        self.assertTrue(res_tie["tied_best_noncognate"])

    def test_per_read_verifier_catches_flag_swap_and_duplicate_ids(self):
        from verify_replication import check_per_read_audit_table

        # Build a valid 2-row synthetic audit table & matching JSON:
        #   r1: cog=1, bt2=1, cog_sim=0.90, all_sim=0.90  (both on)
        #   r2: cog=0, bt2=0, cog_sim=-1.0, all_sim=0.20  (both off)
        header = (
            "read_id\tassigned_amplicon_1based\tis_cognate\tbt2_both_mapped\tbt2_concordant\t"
            "both_mates_ge60\tcog_edlib_sim\tall_edlib_sim\tbest_target_amplicon_1based\t"
            "best_target_strand\tstrictly_better_noncognate\ttied_best_noncognate\n"
        )
        valid_rows = (
            "@r1\t1\t1\t1\t1\t1\t0.900000\t0.900000\t1\t+\t0\t0\n"
            "@r2\t0\t0\t0\t0\t1\t-1.000000\t0.200000\t2\t+\t0\t0\n"
        )
        # Adversarial flag-swap: preserves bt2 marginal (sum=1) and edlib marginal (sum=1),
        # but changes confusion matrix from (1,1,0,0) to (0,0,1,1)!
        swapped_rows = (
            "@r1\t1\t1\t0\t0\t1\t0.900000\t0.900000\t1\t+\t0\t0\n"
            "@r2\t0\t0\t1\t1\t1\t-1.000000\t0.200000\t2\t+\t0\t0\n"
        )
        dup_id_rows = (
            "@r1\t1\t1\t1\t1\t1\t0.900000\t0.900000\t1\t+\t0\t0\n"
            "@r1\t0\t0\t0\t0\t1\t-1.000000\t0.200000\t2\t+\t0\t0\n"
        )
        cm_11 = {"both_ontarget": 1, "both_offtarget": 1, "a_on_b_off": 0, "a_off_b_on": 0}
        audit_json = {
            "all_384_target_concordance_75": {
                **cm_11, "sample_size": 2, "cognate_total": 1,
                "cognate_cross_target_best": 0, "cognate_tied_noncognate_best": 0,
            },
            "all_384_target_concordance_80": {
                **cm_11, "sample_size": 2, "cognate_total": 1,
                "cognate_cross_target_best": 0, "cognate_tied_noncognate_best": 0,
            },
            "matched_pair_bowtie2_audit": {
                "pairs_ge60_both_mates": 2,
                "bt2_both_mapped_any_count": 1,
                "bt2_concordant_count": 1,
                "bt2_both_equals_concordant": True,
                "bt2_ge60_both_mapped_count": 1,
                "vs_cognate_edlib75": dict(cm_11),
                "vs_cognate_edlib80": dict(cm_11),
                "vs_all_target_edlib75": dict(cm_11),
                "vs_all_target_edlib80": dict(cm_11),
            },
        }

        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "audit.tsv.gz"
            with gzip.open(p, "wt", encoding="utf-8", newline="\n") as f:
                f.write(header + valid_rows)
            self.assertTrue(check_per_read_audit_table(p, audit_json))

            with gzip.open(p, "wt", encoding="utf-8", newline="\n") as f:
                f.write(header + swapped_rows)
            self.assertFalse(check_per_read_audit_table(p, audit_json))

            with gzip.open(p, "wt", encoding="utf-8", newline="\n") as f:
                f.write(header + dup_id_rows)
            self.assertFalse(check_per_read_audit_table(p, audit_json))


if __name__ == "__main__":
    unittest.main(verbosity=2)
