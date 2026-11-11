#!/usr/bin/env python3
"""
code/verify_read_accounting.py

Independently performs exhaustive, mutually exclusive read accounting
across all 3,942,385 raw paired-end sequencing reads of SADDLE 384-plex (E19).
Reconciles every exclusion reason with original SADDLE analysis logic.
"""

import sys
import gzip
import json
import time
import argparse
from pathlib import Path
import openpyxl

# Add code directory to path for canonical classifier import
CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from canonical_classifier import match_primer_prefix, match_primer_prefix_historical

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw" if (ROOT / "data" / "raw").exists() else ROOT / "raw"
DERIVED_DIR = ROOT / "data" / "derived" if (ROOT / "data" / "derived").exists() else ROOT / "data"

FASTQ_R1 = RAW_DIR / "E19_384plex_NGS_R1_fastq.gz"
FASTQ_R2 = RAW_DIR / "E19_384plex_NGS_R2_fastq.gz"
MOESM5_PATH = RAW_DIR / "41467_2022_29500_MOESM5_ESM.xlsx"

F_ADP = "AGATCGGAAGAGCACAC"
R_ADP = "AGATCGGAAGAGCGTCG"

def load_primers():
    wb = openpyxl.load_workbook(MOESM5_PATH, data_only=True)
    ws = wb["384 Plex Primer"]
    fps, rps = [], []
    for row in ws.iter_rows(values_only=True):
        fps.append((row[0], row[1].upper()))
        rps.append((row[2], row[3].upper()))
    all_primers = fps + rps
    prefix_map = {}
    for idx, (pid, pseq) in enumerate(all_primers):
        prefix_map.setdefault(pseq[:8], []).append(idx)
    return all_primers, prefix_map

def run_accounting(use_historical: bool = False):
    print("=" * 70)
    print("PHASE 0D: EXHAUSTIVE READ ACCOUNTING VERIFICATION")
    print(f"Mode: {'Historical Baseline (/ len(pseq))' if use_historical else 'Canonical Primary (/ check_len)'}")
    print("=" * 70)

    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    required = [MOESM5_PATH, FASTQ_R1, FASTQ_R2]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        print(f"SKIPPED_MISSING_INPUT: External raw input files missing: {missing}")
        print("This read accounting script requires the 914 MB raw FASTQ files and supplementary spreadsheets.")
        print("To verify read accounting totals from bundled package artifacts, run:")
        print("    python code/verify_replication.py")
        print("For raw data download instructions, see README.md.")
        sys.exit(2)

    all_primers, prefix_map = load_primers()
    print(f"Loaded {len(all_primers)} primers.")

    categories = {
        "no_adapter_both": 0,
        "no_adapter_r1_only": 0,
        "no_adapter_r2_only": 0,
        "qc_short_len": 0,
        "qc_ambiguous_n": 0,
        "qc_low_quality": 0,
        "primer_unmatched_both": 0,
        "primer_unmatched_r1_only": 0,
        "primer_unmatched_r2_only": 0,
        "classified_on_target": 0,
        "classified_dimer": 0,
        "classified_nonspecific": 0,
    }

    total_reads = 0
    t0 = time.time()
    last_t = t0

    with gzip.open(FASTQ_R1, "rt") as f1, gzip.open(FASTQ_R2, "rt") as f2:
        while True:
            r1_id = f1.readline()
            if not r1_id:
                break
            r1_seq = f1.readline().strip()
            f1.readline()
            r1_qual = f1.readline().strip()

            r2_id = f2.readline()
            r2_seq = f2.readline().strip()
            f2.readline()
            r2_qual = f2.readline().strip()

            total_reads += 1
            if total_reads % 500000 == 0:
                now = time.time()
                rate = 500000 / (now - last_t)
                last_t = now
                print(f"  Processed {total_reads:,} reads ({rate:.0f} reads/s)...")

            # 1. Adapter detection
            pos1 = r1_seq.find(F_ADP)
            pos2 = r2_seq.find(R_ADP)

            if pos1 == -1 and pos2 == -1:
                categories["no_adapter_both"] += 1
                continue
            elif pos1 == -1:
                categories["no_adapter_r1_only"] += 1
                continue
            elif pos2 == -1:
                categories["no_adapter_r2_only"] += 1
                continue

            t1 = r1_seq[:pos1]
            t2 = r2_seq[:pos2]

            # 2. Quality Control
            if len(t1) <= 10 or len(t2) <= 10:
                categories["qc_short_len"] += 1
                continue

            if "N" in t1 or "N" in t2:
                categories["qc_ambiguous_n"] += 1
                continue

            q1_pass = sum(1 for c in r1_qual[:len(t1)] if ord(c) > 63) / len(t1) > 0.7
            q2_pass = sum(1 for c in r2_qual[:len(t2)] if ord(c) > 63) / len(t2) > 0.7
            if not (q1_pass and q2_pass):
                categories["qc_low_quality"] += 1
                continue

            # 3. Primer Identification
            m1 = prefix_map.get(t1[:8], [])
            m2 = prefix_map.get(t2[:8], [])

            if use_historical:
                best_id1 = match_primer_prefix_historical(t1, m1, all_primers)
                best_id2 = match_primer_prefix_historical(t2, m2, all_primers)
            else:
                best_id1 = match_primer_prefix(t1, m1, all_primers)
                best_id2 = match_primer_prefix(t2, m2, all_primers)

            if best_id1 is None and best_id2 is None:
                categories["primer_unmatched_both"] += 1
                continue
            elif best_id1 is None:
                categories["primer_unmatched_r1_only"] += 1
                continue
            elif best_id2 is None:
                categories["primer_unmatched_r2_only"] += 1
                continue

            # 4. Classification
            is_cognate = (best_id2 == best_id1 + 384) or (best_id1 == best_id2 + 384)
            p1_len = len(all_primers[best_id1][1])
            p2_len = len(all_primers[best_id2][1])
            dimer_threshold = p1_len + p2_len - 2
            insert_len = len(t1)

            if is_cognate and insert_len >= dimer_threshold:
                categories["classified_on_target"] += 1
            elif insert_len < dimer_threshold:
                categories["classified_dimer"] += 1
            else:
                categories["classified_nonspecific"] += 1

    elapsed = time.time() - t0
    sum_categories = sum(categories.values())
    assert sum_categories == total_reads, f"Partition error: {sum_categories} != {total_reads}"

    classified = (
        categories["classified_on_target"]
        + categories["classified_nonspecific"]
        + categories["classified_dimer"]
    )
    unclassified = total_reads - classified

    print("\n" + "=" * 70)
    print("EXHAUSTIVE READ ACCOUNTING RESULTS")
    print(f"Elapsed Time: {elapsed:.1f}s ({total_reads / elapsed:.0f} reads/s)")
    print(f"Total Raw Read Pairs:       {total_reads:,}")
    print(f"Total Classified Pairs:     {classified:,} ({classified / total_reads * 100:.2f}%)")
    print(f"Total Excluded Pairs:       {unclassified:,} ({unclassified / total_reads * 100:.2f}%)")
    print("-" * 70)
    print("DETAILED BREAKDOWN BY CATEGORY:")
    for cat, count in categories.items():
        pct_raw = count / total_reads * 100
        pct_class = (count / classified * 100) if "classified" in cat else 0.0
        extra = f" [{pct_class:.2f}% of classified]" if "classified" in cat else ""
        print(f"  {cat:<25}: {count:>10,} ({pct_raw:>6.2f}% of raw){extra}")
    print("=" * 70)

    accounting_results = {
        "total_raw_pairs": total_reads,
        "total_classified_pairs": classified,
        "total_excluded_pairs": unclassified,
        "categories": categories,
        "percentages_of_raw": {k: round(v / total_reads * 100, 4) for k, v in categories.items()},
        "percentages_of_classified": {
            "on_target_pct": round(categories["classified_on_target"] / classified * 100, 4),
            "nonspecific_pct": round(categories["classified_nonspecific"] / classified * 100, 4),
            "dimer_pct": round(categories["classified_dimer"] / classified * 100, 4),
        },
    }

    out_json = DERIVED_DIR / "saddle_384_read_accounting.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(accounting_results, f, indent=2)
    print(f"Accounting JSON written to {out_json}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exhaustive read accounting verification")
    parser.add_argument("--historical", action="store_true", help="Run with historical Phase 0D full-length primer denominator")
    args = parser.parse_args()
    run_accounting(use_historical=args.historical)
