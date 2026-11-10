#!/usr/bin/env python3
"""
code/reproduce_from_fastq.py

End-to-End FASTQ-to-Classification Replication Pipeline for:
"[Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing
      Validation of Multiplex PCR Read Partitions"

Author: Rashid Kadayil (ORCID: 0009-0009-6398-4557)
Original Publication: Xie et al., Nature Communications 13, 1881 (2022), Figure 3f
                      DOI: 10.1038/s41467-022-29500-4 | PMID: 35410464

This script processes raw paired-end sequencing reads directly from FASTQ files
(Illumina MiniSeq 2x150 bp, Sample E19, 3,942,385 reads) to replicate the published
read classification breakdown without relying on precomputed tables.

Usage:
    python code/reproduce_from_fastq.py --full       # Complete run across all 3.94M reads (~2 min)
    python code/reproduce_from_fastq.py --sample N   # Fast live test on first N reads (default: 50,000)
"""

import sys
import gzip
import time
import argparse
from pathlib import Path

F_ADP = "AGATCGGAAGAGCACAC"
R_ADP = "AGATCGGAAGAGCGTCG"

def revcomp(s: str) -> str:
    tr = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return s.translate(tr)[::-1]

def find_raw_dir(custom_path=None) -> Path:
    if custom_path:
        return Path(custom_path)
    root = Path(__file__).resolve().parent.parent
    candidates = [
        root / "data" / "raw",
        root.parent / "data" / "raw",
        Path("data/raw"),
        Path("../data/raw"),
        root / "raw",
    ]
    for c in candidates:
        if c.exists() and (c / "41467_2022_29500_MOESM5_ESM.xlsx").exists():
            return c
    return candidates[0]

def load_resources(raw_dir: Path):
    moesm5_path = raw_dir / "41467_2022_29500_MOESM5_ESM.xlsx"
    if not moesm5_path.exists():
        return None, None, None

    try:
        import openpyxl
    except ImportError:
        print("ERROR: openpyxl is required to read primer definitions from MOESM5.", file=sys.stderr)
        sys.exit(1)

    wb = openpyxl.load_workbook(moesm5_path, data_only=True)
    ws_primers = wb["384 Plex Primer"]
    fps, rps = [], []
    for row in ws_primers.iter_rows(values_only=True):
        fps.append((row[0], row[1].upper()))
        rps.append((row[2], row[3].upper()))
    all_primers = fps + rps

    prefix_map = {}
    for idx, (pid, pseq) in enumerate(all_primers):
        prefix_map.setdefault(pseq[:8], []).append(idx)

    # Reference amplicons from 384 Plex Info
    ref_amplicons = {}
    if "384 Plex Info" in wb.sheetnames:
        ws_info = wb["384 Plex Info"]
        for i in range(2, ws_info.max_row + 1):
            fp_id = ws_info[i][0].value
            fp_seq = ws_info[i][1].value.lower()
            rp_id = ws_info[i][2].value
            rp_seq = ws_info[i][3].value.lower()
            full_seq = ws_info[i][9].value.lower() if ws_info[i][9].value else ""
            rc_rp = revcomp(rp_seq)
            pos_fp = full_seq.find(fp_seq)
            pos_rc_rp = full_seq.find(rc_rp)
            if pos_fp != -1 and pos_rc_rp != -1:
                ref_amplicons[i - 2] = full_seq[pos_fp : pos_rc_rp + len(rc_rp)]

    return all_primers, prefix_map, ref_amplicons

def run_pipeline(raw_dir: Path, max_reads: int = None):
    r1_path = raw_dir / "E19_384plex_NGS_R1_fastq.gz"
    r2_path = raw_dir / "E19_384plex_NGS_R2_fastq.gz"
    moesm5_path = raw_dir / "41467_2022_29500_MOESM5_ESM.xlsx"

    required = [r1_path, r2_path, moesm5_path]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        print("=" * 78)
        print("SKIPPED_MISSING_INPUT: External raw FASTQ inputs missing in:")
        print(f"  {raw_dir}")
        print(f"  Missing files: {missing}")
        print("\nTo execute this live FASTQ pipeline, download the raw data from Figshare:")
        print("  DOI: 10.6084/m9.figshare.16944154.v3 (CC BY 4.0)")
        print("  Files required: E19_384plex_NGS_R1_fastq.gz, E19_384plex_NGS_R2_fastq.gz")
        print("And supplementary table MOESM5 from Nature Communications 13, 1881 (2022):")
        print("  File required:  41467_2022_29500_MOESM5_ESM.xlsx")
        print("=" * 78)
        sys.exit(2)

    print("=" * 78)
    print("SADDLE REPLICATION: DIRECT FASTQ-TO-CLASSIFICATION PIPELINE")
    print(f"Raw Input Directory: {raw_dir}")
    print(f"R1: {r1_path.name} ({r1_path.stat().st_size / (1024*1024):.1f} MB)")
    print(f"R2: {r2_path.name} ({r2_path.stat().st_size / (1024*1024):.1f} MB)")
    mode_str = f"FULL DATASET (all 3,942,385 reads)" if max_reads is None else f"SAMPLE BATCH ({max_reads:,} reads)"
    print(f"Execution Mode:      {mode_str}")
    print("=" * 78)

    all_primers, prefix_map, ref_amplicons = load_resources(raw_dir)
    print(f"Loaded {len(all_primers)} primers and {len(ref_amplicons)} reference amplicons.")

    categories = {
        "no_adapter": 0,
        "qc_failed": 0,
        "primer_unmatched": 0,
        "classified_on_target": 0,
        "classified_nonspecific": 0,
        "classified_dimer": 0,
    }

    t0 = time.time()
    last_t = t0
    total_reads = 0

    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
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
                rate = 500000 / (now - last_t) if (now - last_t) > 0 else 0
                last_t = now
                print(f"  Processed {total_reads:,} reads ({rate:.0f} reads/s)...")
                sys.stdout.flush()

            # 1. Adapter detection
            pos1 = r1_seq.find(F_ADP)
            pos2 = r2_seq.find(R_ADP)
            if pos1 == -1 or pos2 == -1:
                categories["no_adapter"] += 1
                if max_reads and total_reads >= max_reads:
                    break
                continue

            t1 = r1_seq[:pos1]
            q1 = r1_qual[:pos1]
            t2 = r2_seq[:pos2]
            q2 = r2_qual[:pos2]

            # 2. Quality filtering (len > 10, no 'N', Q30 > 0.70)
            if len(t1) <= 10 or len(t2) <= 10:
                categories["qc_failed"] += 1
                if max_reads and total_reads >= max_reads:
                    break
                continue
            if "N" in t1 or "N" in t2:
                categories["qc_failed"] += 1
                if max_reads and total_reads >= max_reads:
                    break
                continue

            q1_pass = (sum(1 for c in q1 if ord(c) > 63) / len(q1)) > 0.7
            q2_pass = (sum(1 for c in q2 if ord(c) > 63) / len(q2)) > 0.7
            if not (q1_pass and q2_pass):
                categories["qc_failed"] += 1
                if max_reads and total_reads >= max_reads:
                    break
                continue

            # 3. Primer identification (5' prefix matching)
            best1 = None
            for pid in prefix_map.get(t1[:8], []):
                pseq = all_primers[pid][1]
                check_len = min(len(pseq), len(t1))
                if check_len >= 8 and sum(1 for a, b in zip(t1[:check_len], pseq[:check_len]) if a == b) / check_len > 0.6:
                    best1 = pid
                    break

            best2 = None
            for pid in prefix_map.get(t2[:8], []):
                pseq = all_primers[pid][1]
                check_len = min(len(pseq), len(t2))
                if check_len >= 8 and sum(1 for a, b in zip(t2[:check_len], pseq[:check_len]) if a == b) / check_len > 0.6:
                    best2 = pid
                    break

            if best1 is None or best2 is None:
                categories["primer_unmatched"] += 1
                if max_reads and total_reads >= max_reads:
                    break
                continue

            # 4. Classification
            is_cognate = (best2 == best1 + 384) or (best1 == best2 + 384)
            p1_len = len(all_primers[best1][1])
            p2_len = len(all_primers[best2][1])
            dimer_threshold = p1_len + p2_len - 2
            insert_len = len(t1)

            if is_cognate and insert_len >= dimer_threshold:
                categories["classified_on_target"] += 1
            elif insert_len < dimer_threshold:
                categories["classified_dimer"] += 1
            else:
                categories["classified_nonspecific"] += 1

            if max_reads and total_reads >= max_reads:
                break

    elapsed = time.time() - t0
    classified = categories["classified_on_target"] + categories["classified_nonspecific"] + categories["classified_dimer"]
    excluded = categories["no_adapter"] + categories["qc_failed"] + categories["primer_unmatched"]

    print("-" * 78)
    print("FASTQ PROCESSING RESULTS:")
    print(f"Elapsed Time:               {elapsed:.2f} s ({total_reads / elapsed:.0f} reads/s)")
    print(f"Total Processed Reads:      {total_reads:,}")
    print(f"Total Classified Reads:     {classified:,} ({classified/total_reads*100:.2f}%)")
    print(f"Total Excluded Reads:       {excluded:,} ({excluded/total_reads*100:.2f}%)")
    print(f"  - Missing Adapter:        {categories['no_adapter']:,} ({categories['no_adapter']/total_reads*100:.2f}%)")
    print(f"  - Quality / Length Filter:{categories['qc_failed']:,} ({categories['qc_failed']/total_reads*100:.2f}%)")
    print(f"  - Primer Unmatched:       {categories['primer_unmatched']:,} ({categories['primer_unmatched']/total_reads*100:.2f}%)")

    print("\nCLASSIFIED READ BREAKDOWN:")
    ontarget_pct = categories["classified_on_target"] / classified * 100 if classified > 0 else 0
    nonspec_pct = categories["classified_nonspecific"] / classified * 100 if classified > 0 else 0
    dimer_pct = categories["classified_dimer"] / classified * 100 if classified > 0 else 0

    print(f"  On-Target Amplicons:  {categories['classified_on_target']:>10,} ({ontarget_pct:>6.2f}%) | Published Fig 3f: ~43%")
    print(f"  Nonspecific Products: {categories['classified_nonspecific']:>10,} ({nonspec_pct:>6.2f}%) | Published Fig 3f: ~56%")
    print(f"  Primer-Dimers:        {categories['classified_dimer']:>10,} ({dimer_pct:>6.2f}%) | Published Fig 3f: ~1%")

    if max_reads is None:
        # Full validation checks
        assert total_reads == 3942385, f"Unexpected total reads: {total_reads}"
        assert classified == 2860658, f"Unexpected classified reads: {classified}"
        assert categories["classified_on_target"] == 1254647
        assert categories["classified_nonspecific"] == 1581139
        assert categories["classified_dimer"] == 24872

        print("\n" + "=" * 78)
        print("FULL FASTQ REPLICATION CONFIRMED:")
        print("  All 3,942,385 reads independently processed and verified.")
        print(f"  On-target:   {ontarget_pct:.2f}% (published ~43%, diff = +{ontarget_pct - 43:.2f}%) [PASS within 1%]")
        print(f"  Nonspecific: {nonspec_pct:.2f}% (published ~56%, diff = {nonspec_pct - 56:.2f}%) [PASS within 1%]")
        print(f"  Dimer:       {dimer_pct:.2f}%  (published ~1%,  diff = {dimer_pct - 1:.2f}%)  [PASS within 0.2%]")
        print("=" * 78)
    else:
        print("\n" + "=" * 78)
        print(f"SAMPLE FASTQ TEST COMPLETED ({total_reads:,} reads tested live).")
        print("=" * 78)
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reproduce SADDLE 384-plex Read Accounting from FASTQ")
    parser.add_argument("--full", action="store_true", help="Process the full 3.94M read FASTQ files")
    parser.add_argument("--sample", type=int, default=50000, help="Number of reads for sample check (default: 50,000)")
    parser.add_argument("--raw-dir", type=str, default=None, help="Path to raw FASTQ directory")
    args = parser.parse_args()

    raw_directory = find_raw_dir(args.raw_dir)
    max_r = None if args.full else args.sample
    run_pipeline(raw_directory, max_reads=max_r)
