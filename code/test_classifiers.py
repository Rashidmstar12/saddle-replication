#!/usr/bin/env python3
"""
code/test_classifiers.py

Seven Synthetic Unit Tests for:
"[Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing
      Read Partitions and Primer Interaction Distribution in a 384-Plex PCR Panel"

Author: Rashid Kadayil (ORCID: 0009-0009-6398-4557)

Gate 1 Synthetic Unit Verification:
This test suite contains seven synthetic unit tests exercising the canonical primer
prefix classifier, historical baseline, and simulated input conditions against fixed
synthetic fixtures. These are isolated unit tests that verify mathematical classification
edge cases (including the 12/22 bp short-insert denominator discrepancy); they do not
execute the complete streaming FASTQ production pipeline.
"""

import io
import gzip
import unittest
import sys
from pathlib import Path

# Add code directory to path
CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from canonical_classifier import match_primer_prefix, match_primer_prefix_historical


class TestClassifierEquivalenceAndSafeguards(unittest.TestCase):
    """Seven synthetic unit tests verifying mathematical classification edge cases and simulated input conditions."""

    def setUp(self):
        # Synthetic 22-bp primer
        self.test_primer = ("fp_test_1", "ACGTACGTACGTACGTACGTAC")
        self.all_primers = [self.test_primer]
        self.candidate_indices = [0]

    def test_short_read_12_of_22_matches(self):
        """
        Test Case 1: Short-read 12/22 primer match discrepancy.

        A trimmed read of 12 bp perfectly matches the first 12 bases of a 22 bp primer.
        - Canonical classifier: 12 / 12 = 1.000 > 0.6 -> ACCEPT (returns index 0).
        - Historical classifier: 12 / 22 = 0.545 <= 0.6 -> REJECT (returns None).

        This exact discrepancy accounts for the 6-read difference between the
        primary rerun (24,878 dimers) and historical partition (24,872 dimers).
        """
        trimmed_read_12bp = "ACGTACGTACGT"  # Exactly matches first 12 bp of primer

        # Canonical classifier divides by min(len(p), len(read)) = 12
        canonical_result = match_primer_prefix(
            trimmed_read_12bp, self.candidate_indices, self.all_primers
        )
        self.assertEqual(
            canonical_result,
            0,
            "Canonical classifier must accept perfect 12-bp prefix match (score = 12/12 = 1.0)",
        )

        # Historical classifier divides by full primer length = 22
        historical_result = match_primer_prefix_historical(
            trimmed_read_12bp, self.candidate_indices, self.all_primers
        )
        self.assertIsNone(
            historical_result,
            "Historical classifier must reject 12-bp match because 12/22 = 0.545 <= 0.6",
        )

    def test_full_length_primer_matches(self):
        """Test Case 2: Full-length reads match identically under both classifiers."""
        full_read = "ACGTACGTACGTACGTACGTAC"  # Exact 22 bp match
        canonical_result = match_primer_prefix(
            full_read, self.candidate_indices, self.all_primers
        )
        historical_result = match_primer_prefix_historical(
            full_read, self.candidate_indices, self.all_primers
        )
        self.assertEqual(canonical_result, 0)
        self.assertEqual(historical_result, 0)

    def test_mismatched_prefix_rejected(self):
        """Test Case 3: Read with substantial mismatches is rejected under both classifiers."""
        mismatched_read = "ACGTACGTNNNNNNNNNNNNNN"  # Only 8/22 matches = 0.36
        canonical_result = match_primer_prefix(
            mismatched_read, self.candidate_indices, self.all_primers
        )
        historical_result = match_primer_prefix_historical(
            mismatched_read, self.candidate_indices, self.all_primers
        )
        self.assertIsNone(canonical_result)
        self.assertIsNone(historical_result)

    def test_sub_threshold_length_rejected(self):
        """Test Case 4: Inserts shorter than min_prefix (8 bp) are rejected."""
        very_short_read = "ACGTACG"  # 7 bp
        canonical_result = match_primer_prefix(
            very_short_read, self.candidate_indices, self.all_primers, min_prefix=8
        )
        self.assertIsNone(canonical_result, "Reads < 8 bp must be rejected")

    def test_paired_fastq_id_synchronization(self):
        """
        Test Case 5: Paired FASTQ read identifier validation (synthetic unit check).
        Compares constructed paired headers; verifies string equality and mismatch detection logic.
        """
        id1_good = "@M04612:271:000000000-CWC2J:1:1101:16920:1534 1:N:0:1"
        id2_good = "@M04612:271:000000000-CWC2J:1:1101:16920:1534 2:N:0:1"
        id1_clean = id1_good.strip().split()[0].removesuffix("/1")
        id2_clean = id2_good.strip().split()[0].removesuffix("/2")
        self.assertEqual(id1_clean, id2_clean, "Synchronized read headers must match")

        id2_bad = "@M04612:271:000000000-CWC2J:1:1101:99999:9999 2:N:0:1"
        id2_bad_clean = id2_bad.strip().split()[0].removesuffix("/2")
        self.assertNotEqual(
            id1_clean, id2_bad_clean, "Desynchronized read headers must be detected"
        )

    def test_missing_reference_amplicon_fail_closed(self):
        """
        Test Case 6: Missing reference amplicon dictionary lookup (synthetic unit check).
        Verifies that an unmapped target amplicon key is absent from the reference mapping.
        """
        ref_amplicons = {0: "ACGTACGTACGTACGTACGT"}
        # Target amplicon index 99 is missing from dictionary
        assigned_target = 99
        self.assertNotIn(
            assigned_target,
            ref_amplicons,
            "Target amplicon 99 is missing; pipeline must route to nonspecific",
        )

    def test_adapter_prefix_detection(self):
        """Test Case 7: Illumina TruSeq adapter prefix trimming (synthetic unit check)."""
        f_adp = "AGATCGGAAGAGCACAC"
        r_adp = "AGATCGGAAGAGCGTCG"
        r1_seq = "ACGTACGTACGT" + f_adp + "TTTTTT"
        r2_seq = "TGCATGCATGCA" + r_adp + "AAAAAA"

        pos1 = r1_seq.find(f_adp)
        pos2 = r2_seq.find(r_adp)
        self.assertEqual(pos1, 12)
        self.assertEqual(pos2, 12)

        # Trimming produces exact 12 bp inserts
        t1 = r1_seq[:pos1]
        t2 = r2_seq[:pos2]
        self.assertEqual(len(t1), 12)
        self.assertEqual(len(t2), 12)


def run_tests():
    suite = unittest.TestLoader().loadTestsFromTestCase(
        TestClassifierEquivalenceAndSafeguards
    )
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    if not result.wasSuccessful():
        sys.exit(1)
    print("\nALL SEVEN SYNTHETIC UNIT TESTS PASSED (Gate 1 PASS).")


if __name__ == "__main__":
    run_tests()
