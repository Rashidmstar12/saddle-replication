#!/usr/bin/env python3
"""
code/run_random_audits.py

Uniform random-sample audits (seeded reservoir sampling across the entire
3,942,385-read E19 library, eliminating flow-cell tile order bias):
  1. Matched-pair alignment audit on random N=5,000 qualifying Pipeline B
     non-dimer read pairs (seed=20261124):
       Evaluates ALL THREE alignment strategies on the EXACT SAME 5,000 pairs:
         (a) Cognate-anchored edlib (>=75% and >=80%)
         (b) Dual-strand all-384-target edlib (>=75% and >=80%, 768 targets/read)
         (c) Paired-end Bowtie2 against the exact 384-amplicon reference
       Archives the complete per-read evidence table to
       data/matched_pair_audit_5000.tsv.gz (tracking both strictly better AND
       tied non-cognate best matches) and reports exact 2x2 confusion matrices
       and Wilson 95% confidence intervals.
  2. Gapped vs. ungapped indel CIGAR recovery audit on random N=20,000
     qualifying cognate non-dimer reads (seed=20261125).
"""

import argparse
import gzip
import json
import math
import random
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import edlib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import canonical_classifier as cc
from run_saddle_replication import (
    PANELS,
    find_raw_dir,
    load_panel_primers_and_refs,
    require_bowtie2,
)


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1.0 + (z * z) / n
    center = (p + (z * z) / (2.0 * n)) / denom
    half = (z / denom) * math.sqrt((p * (1.0 - p)) / n + (z * z) / (4.0 * n * n))
    return round(max(0.0, center - half) * 100.0, 4), round(min(1.0, center + half) * 100.0, 4)


def parse_cigar_indels(cigar: str) -> Tuple[int, int, int, int]:
    ins = dels = mis = mat = 0
    for count_s, op in re.findall(r"(\d+)([=XID])", cigar):
        c = int(count_s)
        if op == "I":
            ins += c
        elif op == "D":
            dels += c
        elif op == "X":
            mis += c
        elif op == "=":
            mat += c
    return ins, dels, mis, mat


def run_bowtie2_on_sample(
    sample_records: List[Tuple[str, str, str, str, str, str, int]],
    ref_fasta: Path,
    bt2: str,
    bt2_build: str,
    threads: int = 4,
) -> List[Tuple[bool, bool, bool]]:
    n = len(sample_records)
    results: List[Tuple[bool, bool, bool]] = [(False, False, False)] * n
    with tempfile.TemporaryDirectory(prefix="saddle_audit_bt2_") as tmp_str:
        tmp = Path(tmp_str)
        r1_fq = tmp / "sample_R1.fastq"
        r2_fq = tmp / "sample_R2.fastq"
        idx_prefix = tmp / "ref_idx"
        sam_path = tmp / "sample.sam"

        with (
            open(r1_fq, "w", encoding="ascii", newline="\n") as out1,
            open(r2_fq, "w", encoding="ascii", newline="\n") as out2,
        ):
            for idx, (_, t1, tq1, t2, tq2, _, _) in enumerate(sample_records):
                out1.write(f"@s{idx}\n{t1}\n+\n{tq1}\n")
                out2.write(f"@s{idx}\n{t2}\n+\n{tq2}\n")

        subprocess.run(
            [bt2_build, "--quiet", str(ref_fasta), str(idx_prefix)],
            check=True, capture_output=True, text=True,
        )
        subprocess.run(
            [bt2, "-x", str(idx_prefix), "-1", str(r1_fq), "-2", str(r2_fq), "--threads", str(threads), "-S", str(sam_path)],
            check=True, capture_output=True, text=True,
        )

        pending: Dict[str, int] = {}
        with open(sam_path, "r", encoding="ascii") as sf:
            for line in sf:
                if not line or line[0] == "@":
                    continue
                parts = line.split("\t", 3)
                qname = parts[0]
                flag = int(parts[1])
                if flag & 0x900:
                    continue
                if qname not in pending:
                    pending[qname] = flag
                else:
                    f_prev = pending.pop(qname)
                    idx = int(qname[1:])
                    r1_unm = bool(f_prev & 0x4)
                    r2_unm = bool(flag & 0x4)
                    both_mapped = (not r1_unm) and (not r2_unm)
                    concordant = both_mapped and bool(f_prev & 0x2)
                    _, t1, _, t2, _, _, _ = sample_records[idx]
                    ge60 = (len(t1) >= 60) and (len(t2) >= 60)
                    results[idx] = (both_mapped, concordant, ge60)
    return results


def build_binary_cm(pred_a: List[bool], pred_b: List[bool], name_a: str, name_b: str) -> dict:
    assert len(pred_a) == len(pred_b)
    n = len(pred_a)
    both_on = sum(1 for a, b in zip(pred_a, pred_b) if a and b)
    both_off = sum(1 for a, b in zip(pred_a, pred_b) if (not a) and (not b))
    a_on_b_off = sum(1 for a, b in zip(pred_a, pred_b) if a and (not b))
    a_off_b_on = sum(1 for a, b in zip(pred_a, pred_b) if (not a) and b)
    agree = both_on + both_off
    disc = a_on_b_off + a_off_b_on
    assert agree + disc == n
    a_on = both_on + a_on_b_off
    b_on = both_on + a_off_b_on
    return {
        "sample_size": n,
        "method_a": name_a,
        "method_b": name_b,
        "both_ontarget": both_on,
        "both_offtarget": both_off,
        "a_on_b_off": a_on_b_off,
        "a_off_b_on": a_off_b_on,
        "method_a_ontarget_pct": round(a_on / n * 100.0, 4),
        "method_a_ontarget_ci95": list(wilson_ci(a_on, n)),
        "method_b_ontarget_pct": round(b_on / n * 100.0, 4),
        "method_b_ontarget_ci95": list(wilson_ci(b_on, n)),
        "diff_b_minus_a_pp": round((b_on - a_on) / n * 100.0, 4),
        "concordance_pct": round(agree / n * 100.0, 4),
        "concordance_ci95": list(wilson_ci(agree, n)),
        "discordance_pct": round(disc / n * 100.0, 4),
        "discordance_ci95": list(wilson_ci(disc, n)),
    }


def run_audits(
    raw_dir: Path,
    out_data_dir: Path,
    seed: int = 20261124,
    n_concordance: int = 5000,
    n_indel: int = 20000,
    threads: int = 4,
    skip_bowtie2: bool = False,
) -> dict:
    run_bt2, bt2, bt2_build, _ = require_bowtie2(skip_bowtie2)

    base_dir = Path(__file__).resolve().parent.parent
    source_ref_dir = base_dir / "data" / "reference_fastas"
    out_ref_dir = out_data_dir / "reference_fastas"
    out_ref_dir.mkdir(parents=True, exist_ok=True)
    moesm5_path = raw_dir / "41467_2022_29500_MOESM5_ESM.xlsx"
    e19_cfg = [p for p in PANELS if p["panel_id"] == "E19"][0]

    _, primers_lower, prefix_index, amplicons, ref_fasta = load_panel_primers_and_refs(
        e19_cfg, moesm5_path, source_ref_dir, out_ref_dir
    )
    amplicons_rc = [cc.revcomp(a) for a in amplicons]
    plex = 384

    r1_path = raw_dir / e19_cfg["r1"]
    r2_path = raw_dir / e19_cfg["r2"]

    rng_conc = random.Random(seed)
    rng_indel = random.Random(seed + 1)

    res_conc: List[Tuple[str, str, str, str, str, str, int]] = []
    res_indel: List[Tuple[str, str, int]] = []

    seen_nondimer = 0
    seen_cognate = 0
    t0 = time.time()
    print(f"Reservoir-sampling full validated E19 library (seed={seed})...", flush=True)

    for rid, s1, q1, s2, q2 in cc.iter_validated_fastq_pairs(r1_path, r2_path):
        t1, tq1, t2, tq2, adp_st = cc.trim_adapter_pair(
            s1, q1, s2, q2, cc.CORE_ADAPTER_R1, cc.CORE_ADAPTER_R2
        )
        if adp_st != "ok":
            continue
        if cc.qc_check_pair(t1, tq1, t2, tq2) != "ok":
            continue

        t1_low = t1.lower()
        t2_low = t2.lower()
        id1 = cc.match_primer_prefix(t1_low, prefix_index.get(t1_low[:8], []), primers_lower)
        if id1 is None:
            continue
        id2 = cc.match_primer_prefix(t2_low, prefix_index.get(t2_low[:8], []), primers_lower)
        if id2 is None:
            continue
        if cc.is_dimer_insert(len(t1_low), len(primers_lower[id1]), len(primers_lower[id2])):
            continue

        seen_nondimer += 1
        amp_idx = -1
        orient_read = t1_low
        if id1 < plex and id2 == id1 + plex:
            amp_idx = id1
            orient_read = t1_low
        elif id2 < plex and id1 == id2 + plex:
            amp_idx = id2
            orient_read = t2_low

        if len(res_conc) < n_concordance:
            res_conc.append((rid, t1, tq1, t2, tq2, orient_read, amp_idx))
        else:
            j = rng_conc.randint(0, seen_nondimer - 1)
            if j < n_concordance:
                res_conc[j] = (rid, t1, tq1, t2, tq2, orient_read, amp_idx)

        if amp_idx != -1:
            seen_cognate += 1
            if len(res_indel) < n_indel:
                res_indel.append((rid, orient_read, amp_idx))
            else:
                j2 = rng_indel.randint(0, seen_cognate - 1)
                if j2 < n_indel:
                    res_indel[j2] = (rid, orient_read, amp_idx)

    print(
        f"Sampled {len(res_conc)} non-dimer (from {seen_nondimer:,}) and "
        f"{len(res_indel)} cognate (from {seen_cognate:,}) in {time.time() - t0:.1f}s.",
        flush=True,
    )

    # Run Bowtie2 first on the exact same 5,000 sampled pairs
    bt2_res = run_bowtie2_on_sample(res_conc, ref_fasta, bt2, bt2_build, threads=threads) if run_bt2 else [(False, False, False)] * len(res_conc)

    # Audit 1: Per-read evaluation across all 384 amplicons (both strands collapsed per amplicon)
    t1 = time.time()
    cog75_flags: List[bool] = []
    cog80_flags: List[bool] = []
    all75_flags: List[bool] = []
    all80_flags: List[bool] = []
    cognate_total = 0
    cognate_strict_cross_75 = 0
    cognate_strict_cross_80 = 0
    cognate_tied_cross_75 = 0
    cognate_tied_cross_80 = 0

    tsv_lines = [
        "read_id\tassigned_amplicon_1based\tis_cognate\tbt2_both_mapped\tbt2_concordant\t"
        "both_mates_ge60\tcog_edlib_sim\tall_edlib_sim\tbest_target_amplicon_1based\t"
        "best_target_strand\tstrictly_better_noncognate\ttied_best_noncognate\n"
    ]

    for idx_r, (rid, t1_seq, _, t2_seq, _, rseq, amp_idx) in enumerate(res_conc):
        cls = cc.classify_best_target_edlib(rseq, amplicons, amplicons_rc, amp_idx)
        assigned_sim = cls["assigned_fwd_sim"]
        best_sim = cls["best_overall_sim"]
        best_amp = cls["best_amp_1based"]
        best_strand = cls["best_strand"]
        strict_better_noncog = cls["strictly_better_noncognate"]
        tied_noncognate = cls["tied_best_noncognate"]

        c75 = (amp_idx != -1) and (assigned_sim >= 0.75)
        c80 = (amp_idx != -1) and (assigned_sim >= 0.80)
        a75 = best_sim >= 0.75
        a80 = best_sim >= 0.80

        cog75_flags.append(c75)
        cog80_flags.append(c80)
        all75_flags.append(a75)
        all80_flags.append(a80)

        if amp_idx != -1:
            cognate_total += 1
            if a75:
                if strict_better_noncog:
                    cognate_strict_cross_75 += 1
                elif tied_noncognate:
                    cognate_tied_cross_75 += 1
            if a80:
                if strict_better_noncog:
                    cognate_strict_cross_80 += 1
                elif tied_noncognate:
                    cognate_tied_cross_80 += 1

        bt2_any, bt2_conc, ge60 = bt2_res[idx_r]
        if not run_bt2:
            ge60 = (len(t1_seq) >= 60) and (len(t2_seq) >= 60)
        tsv_lines.append(
            f"{rid}\t{amp_idx + 1 if amp_idx != -1 else 0}\t{1 if amp_idx != -1 else 0}\t"
            f"{1 if bt2_any else 0}\t{1 if bt2_conc else 0}\t{1 if ge60 else 0}\t"
            f"{assigned_sim:.6f}\t{best_sim:.6f}\t{best_amp}\t{best_strand}\t"
            f"{1 if strict_better_noncog else 0}\t{1 if tied_noncognate else 0}\n"
        )

        if (idx_r + 1) % 1000 == 0:
            print(f"  Matched-pair edlib audit: {idx_r + 1}/{len(res_conc)} ({time.time() - t1:.1f}s)", flush=True)

    # Save deterministic gzipped per-read evidence table (mtime=0)
    per_read_tsv = out_data_dir / "matched_pair_audit_5000.tsv.gz"
    per_read_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(per_read_tsv, "wb") as raw_f:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_f, mtime=0) as gz:
            gz.write("".join(tsv_lines).encode("utf-8"))
    print(f"Wrote 5,000-pair per-read evidence table -> {per_read_tsv}", flush=True)

    bt2_matched = None
    if run_bt2:
        bt2_any_flags = [x[0] for x in bt2_res]
        bt2_conc_flags = [x[1] for x in bt2_res]
        bt2_ge60_any_flags = [x[0] and x[2] for x in bt2_res]
        n_ge60 = sum(1 for x in bt2_res if x[2])
        n_bt2_any = sum(bt2_any_flags)
        n_bt2_conc = sum(bt2_conc_flags)
        n_bt2_ge60 = sum(bt2_ge60_any_flags)
        bt2_matched = {
            "sample_size": len(res_conc),
            "pairs_ge60_both_mates": n_ge60,
            "pairs_lt60_either_mate": len(res_conc) - n_ge60,
            "bt2_both_mapped_any_count": n_bt2_any,
            "bt2_both_mapped_any_pct": round(n_bt2_any / len(res_conc) * 100.0, 4),
            "bt2_both_mapped_any_ci95": list(wilson_ci(n_bt2_any, len(res_conc))),
            "bt2_concordant_count": n_bt2_conc,
            "bt2_concordant_pct": round(n_bt2_conc / len(res_conc) * 100.0, 4),
            "bt2_concordant_ci95": list(wilson_ci(n_bt2_conc, len(res_conc))),
            "bt2_both_equals_concordant": bool(bt2_any_flags == bt2_conc_flags),
            "bt2_ge60_both_mapped_count": n_bt2_ge60,
            "bt2_ge60_both_mapped_pct": round(n_bt2_ge60 / len(res_conc) * 100.0, 4),
            "bt2_ge60_both_mapped_ci95": list(wilson_ci(n_bt2_ge60, len(res_conc))),
            "vs_cognate_edlib75": build_binary_cm(bt2_any_flags, cog75_flags, "bowtie2_paired", "cognate_edlib75"),
            "vs_cognate_edlib80": build_binary_cm(bt2_any_flags, cog80_flags, "bowtie2_paired", "cognate_edlib80"),
            "vs_all_target_edlib75": build_binary_cm(bt2_any_flags, all75_flags, "bowtie2_paired", "all_target_edlib75"),
            "vs_all_target_edlib80": build_binary_cm(bt2_any_flags, all80_flags, "bowtie2_paired", "all_target_edlib80"),
        }

    conc75 = build_binary_cm(cog75_flags, all75_flags, "cognate_edlib75", "all_target_edlib75")
    conc75["cognate_total"] = cognate_total
    conc75["noncognate_total"] = len(res_conc) - cognate_total
    conc75["cognate_cross_target_best"] = cognate_strict_cross_75
    conc75["cognate_tied_noncognate_best"] = cognate_tied_cross_75

    conc80 = build_binary_cm(cog80_flags, all80_flags, "cognate_edlib80", "all_target_edlib80")
    conc80["cognate_total"] = cognate_total
    conc80["noncognate_total"] = len(res_conc) - cognate_total
    conc80["cognate_cross_target_best"] = cognate_strict_cross_80
    conc80["cognate_tied_noncognate_best"] = cognate_tied_cross_80

    # Audit 2: Indel CIGAR recovery audit on random N=20,000 cognate reads
    t2 = time.time()
    recovered_records = []
    passed_gapped = 0
    passed_ungapped = 0
    for rid, rseq, amp_idx in res_indel:
        ref = amplicons[amp_idx]
        g_sim = cc.edlib_shw_identity(rseq, ref)
        u_sim = cc.ungapped_positional_identity(rseq, ref)
        if g_sim >= 0.75:
            passed_gapped += 1
        if u_sim >= 0.75:
            passed_ungapped += 1
        if g_sim >= 0.75 and u_sim < 0.75:
            aln = edlib.align(rseq, ref, mode="SHW", task="path")
            cigar = aln["cigar"]
            ins, dels, mis, mat = parse_cigar_indels(cigar)
            recovered_records.append({
                "read_id": rid,
                "amplicon_index_1based": amp_idx + 1,
                "read_len": len(rseq),
                "gapped_sim": round(g_sim, 4),
                "ungapped_sim": round(u_sim, 4),
                "insertions_bp": ins,
                "deletions_bp": dels,
                "total_indel_bp": ins + dels,
                "substitutions_bp": mis,
                "cigar": cigar,
            })

    n_rec = len(recovered_records)
    indel_totals = np.array([r["total_indel_bp"] for r in recovered_records], dtype=np.int64)
    n_with_gaps = int((indel_totals > 0).sum()) if n_rec else 0
    n_le2 = int(((indel_totals > 0) & (indel_totals <= 2)).sum()) if n_rec else 0
    n_gt2 = int((indel_totals > 2).sum()) if n_rec else 0

    indel_summary = {
        "sampling_method": "uniform_reservoir_over_all_E19_cognate_nondimer_reads",
        "seed": seed + 1,
        "population_cognate_nondimer_reads": seen_cognate,
        "sample_size": len(res_indel),
        "passed_gapped_75": passed_gapped,
        "passed_ungapped_75": passed_ungapped,
        "recovered_by_gapped_count": n_rec,
        "recovered_pct_of_cognate_sample": round(n_rec / len(res_indel) * 100.0, 4),
        "recovered_pct_ci95": list(wilson_ci(n_rec, len(res_indel))),
        "reads_with_gaps_count": n_with_gaps,
        "reads_with_gaps_pct": round(n_with_gaps / n_rec * 100.0, 4) if n_rec else 0.0,
        "substitution_only_count": n_rec - n_with_gaps,
        "small_indels_le2bp_count": n_le2,
        "small_indels_le2bp_pct": round(n_le2 / n_rec * 100.0, 4) if n_rec else 0.0,
        "larger_or_multi_indels_gt2bp_count": n_gt2,
        "larger_or_multi_indels_gt2bp_pct": round(n_gt2 / n_rec * 100.0, 4) if n_rec else 0.0,
        "median_indel_bp": float(np.median(indel_totals)) if n_rec else 0.0,
        "iqr_indel_bp": [
            float(np.percentile(indel_totals, 25)),
            float(np.percentile(indel_totals, 75)),
        ] if n_rec else [0.0, 0.0],
        "elapsed_seconds": round(time.time() - t2, 2),
    }

    audit_out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "seed": seed,
        "population_nondimer_reads": seen_nondimer,
        "population_cognate_nondimer_reads": seen_cognate,
        "per_read_evidence_table": "matched_pair_audit_5000.tsv.gz",
        "all_384_target_concordance_75": conc75,
        "all_384_target_concordance_80": conc80,
        "matched_pair_bowtie2_audit": bt2_matched,
        "indel_recovery_audit": indel_summary,
        "recovered_cigar_records": recovered_records,
    }

    audit_path = out_data_dir / "random_audit_results.json"
    cc.write_utf8_lf(audit_path, json.dumps(audit_out, indent=2) + "\n")
    print(f"Wrote random audits -> {audit_path}", flush=True)

    canon_path = out_data_dir / "canonical_results.json"
    if canon_path.exists():
        canon = json.loads(canon_path.read_text(encoding="utf-8"))
        canon["random_audits_e19"] = {
            "seed": seed,
            "population_nondimer_reads": seen_nondimer,
            "population_cognate_nondimer_reads": seen_cognate,
            "per_read_evidence_table": "matched_pair_audit_5000.tsv.gz",
            "all_384_target_concordance_75": conc75,
            "all_384_target_concordance_80": conc80,
            "matched_pair_bowtie2_audit": bt2_matched,
            "indel_recovery_audit": indel_summary,
        }
        cc.write_utf8_lf(canon_path, json.dumps(canon, indent=2) + "\n")
        print(f"Merged matched-pair random audit summary into -> {canon_path}", flush=True)

    return audit_out


def main() -> None:
    parser = argparse.ArgumentParser(description="Run seeded uniform matched-pair random audits on E19.")
    parser.add_argument("--raw-dir", type=str, default=None)
    parser.add_argument("--out-dir", type=str, default=None)
    parser.add_argument("--seed", type=int, default=20261124)
    parser.add_argument("--n-concordance", type=int, default=5000)
    parser.add_argument("--n-indel", type=int, default=20000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip-bowtie2", action="store_true")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent.parent
    raw_dir = find_raw_dir(base_dir, args.raw_dir)
    out_data_dir = Path(args.out_dir).resolve() if args.out_dir else (base_dir / "data")
    run_audits(
        raw_dir, out_data_dir, seed=args.seed,
        n_concordance=args.n_concordance, n_indel=args.n_indel,
        threads=args.threads, skip_bowtie2=args.skip_bowtie2,
    )


if __name__ == "__main__":
    main()
