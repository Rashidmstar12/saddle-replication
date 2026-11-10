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

It performs a dual methodological evaluation:
1. Method 1 (Primer-Pairing Heuristic): Evaluates cognate primer pairing (fp_i x rp_i)
   and insert length >= L1 + L2 - 2 bp.
2. Method 2 (Cognate-Anchored Gapped Sequence Alignment via edlib): Executes Myers' exact
   bit-vector semi-global alignment (mode='SHW') against the 384 intended reference
   amplicons from Supplementary Data MOESM5.

Strict Fail-Closed Validations:
- Requires edlib for gapped alignment; fails with exit code 1 if missing.
- Enforces paired-end FASTQ read identifier matching; raises ValueError on mismatch.
- Rejects missing reference amplicons as non-target; requires all 384 reference amplicons.
- Computes PASS/MATCH status dynamically based on quantitative tolerance.

Usage:
    python code/reproduce_from_fastq.py --full       # Complete run across all 3.94M reads (~2 min)
    python code/reproduce_from_fastq.py --sample N   # Fast live test on first N reads (default: 50,000, ~1 s)
"""

import sys
import gzip
import time
import argparse
import hashlib
import platform
import datetime
from pathlib import Path

# Defect B Fix: Require edlib explicitly with no silent fallback
try:
    import edlib
except ImportError:
    print("ERROR: 'edlib' is required for genuine gapped sequence alignment (Method 2).", file=sys.stderr)
    print("Please install via: pip install edlib", file=sys.stderr)
    sys.exit(1)

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

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
        print("Install via: pip install openpyxl", file=sys.stderr)
        sys.exit(1)

    wb = openpyxl.load_workbook(moesm5_path, data_only=True)
    ws_primers = wb["384 Plex Primer"]
    fps, rps = [], []
    for row in ws_primers.iter_rows(values_only=True):
        fps.append((row[0], row[1].upper()))
        rps.append((row[2], row[3].upper()))
    all_primers = fps + rps

    if len(all_primers) != 768:
        raise ValueError(f"Expected 768 primers in MOESM5 '384 Plex Primer', found {len(all_primers)}")

    prefix_map = {}
    for idx, (pid, pseq) in enumerate(all_primers):
        prefix_map.setdefault(pseq[:8], []).append(idx)

    # Reference amplicons from 384 Plex Info
    ref_amplicons = {}
    if "384 Plex Info" in wb.sheetnames:
        ws_info = wb["384 Plex Info"]
        for i in range(2, ws_info.max_row + 1):
            fp_seq = ws_info[i][1].value.lower() if ws_info[i][1].value else ""
            rp_seq = ws_info[i][3].value.lower() if ws_info[i][3].value else ""
            full_seq = ws_info[i][9].value.lower() if ws_info[i][9].value else ""
            rc_rp = revcomp(rp_seq)
            pos_fp = full_seq.find(fp_seq)
            pos_rc_rp = full_seq.find(rc_rp)
            if pos_fp != -1 and pos_rc_rp != -1:
                ref_amplicons[i - 2] = full_seq[pos_fp : pos_rc_rp + len(rc_rp)]

    if len(ref_amplicons) != 384:
        raise ValueError(f"Expected 384 reference amplicons in MOESM5 '384 Plex Info', found {len(ref_amplicons)}")

    return all_primers, prefix_map, ref_amplicons

def run_pipeline(raw_dir: Path, max_reads: int = None):
    r1_path = raw_dir / "E19_384plex_NGS_R1_fastq.gz"
    r2_path = raw_dir / "E19_384plex_NGS_R2_fastq.gz"
    moesm5_path = raw_dir / "41467_2022_29500_MOESM5_ESM.xlsx"

    required = [r1_path, r2_path, moesm5_path]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        print("\n" + "=" * 78, file=sys.stderr)
        print("ERROR: Missing required raw data input files:", file=sys.stderr)
        print(f"  Directory checked: {raw_dir}", file=sys.stderr)
        print(f"  Missing files:     {missing}", file=sys.stderr)
        print("\nTo run this script, download the raw sequencing files from Figshare:", file=sys.stderr)
        print("  DOI: 10.6084/m9.figshare.16944154.v3 (CC BY 4.0)", file=sys.stderr)
        print("  Files: E19_384plex_NGS_R1_fastq.gz, E19_384plex_NGS_R2_fastq.gz", file=sys.stderr)
        print("And supplementary table MOESM5 from Nature Communications 13, 1881 (2022):", file=sys.stderr)
        print("  File:  41467_2022_29500_MOESM5_ESM.xlsx", file=sys.stderr)
        print("=" * 78, file=sys.stderr)
        sys.exit(2)

    import openpyxl
    all_primers, prefix_map, ref_amplicons = load_resources(raw_dir)
    if all_primers is None:
        print("ERROR: Failed to load primer sequences.", file=sys.stderr)
        sys.exit(2)

    mode_str = f"Streaming first {max_reads:,} reads" if max_reads else "Full run (all 3,942,385 reads)"
    print("=" * 78)
    print("RESCIENCE C SADDLE FASTQ REPLICATION PIPELINE - EXECUTION AUDIT")
    print(f"Timestamp:       {datetime.datetime.now(datetime.timezone.utc).isoformat()}")
    print(f"Platform:        {platform.platform()}")
    print(f"Python:          {sys.version.split()[0]} ({sys.executable})")
    print(f"openpyxl:        {openpyxl.__version__}")
    print(f"edlib:           Myers' exact bit-vector semi-global alignment (SHW)")
    print(f"Dataset:         Sample E19 (384-plex, Illumina MiniSeq 2x150 bp)")
    print(f"Mode:            {mode_str}")
    print("Input File Checksums (SHA-256):")
    print(f"  {moesm5_path.name}: {sha256_file(moesm5_path)}")
    print(f"  {r1_path.name}: {sha256_file(r1_path)}")
    print(f"  {r2_path.name}: {sha256_file(r2_path)}")
    print("=" * 78)

    t0 = time.time()
    total_reads = 0
    categories = {
        "no_adapter": 0,
        "qc_failed": 0,
        "primer_unmatched": 0,
        "dimers": 0,
        "method1_ontarget": 0,
        "method1_nonspecific": 0,
        "method2_ontarget_75": 0,
        "method2_nonspecific_75": 0,
        "method2_ontarget_80": 0,
        "method2_nonspecific_80": 0,
    }

    last_time = t0
    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
        while True:
            r1_id = f1.readline()
            if not r1_id:
                break
            r1_seq = f1.readline().strip()
            f1.readline()
            r1_qual = f1.readline().strip()

            r2_id = f2.readline()
            if not r2_id:
                raise ValueError(f"FASTQ pair truncation: R1 has more reads than R2 at read {total_reads+1}")
            r2_seq = f2.readline().strip()
            f2.readline()
            r2_qual = f2.readline().strip()

            total_reads += 1

            # Validate FASTQ paired-end read identifier synchronization
            id1 = r1_id.strip().split()[0].removesuffix("/1")
            id2 = r2_id.strip().split()[0].removesuffix("/2")
            if id1 != id2:
                raise ValueError(f"FASTQ pair desynchronization detected at read {total_reads}: R1 '{id1}' != R2 '{id2}'")

            if total_reads % 500000 == 0:
                cur_t = time.time()
                rate = 500000 / (cur_t - last_time) if cur_t > last_time else 0
                last_time = cur_t
                print(f"  Processed {total_reads:,} reads ({rate:.0f} reads/s)...")

            # 1. Adapter identification and trimming
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
            if len(t1) <= 10 or len(t2) <= 10 or "N" in t1 or "N" in t2:
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

            if insert_len < dimer_threshold:
                categories["dimers"] += 1
            elif is_cognate:
                # Method 1: Primer Pairing
                categories["method1_ontarget"] += 1

                # Method 2: Cognate-Anchored Gapped Sequence Alignment via edlib
                amp_idx = best1 if best1 < 384 else best2
                exp_amp = ref_amplicons.get(amp_idx, "")
                if best1 >= 384 and exp_amp:
                    exp_amp = revcomp(exp_amp)

                # Defect A Fix: Fail-closed on missing reference amplicon sequence
                if exp_amp:
                    res = edlib.align(t1.lower(), exp_amp, mode="SHW", task="distance")
                    # Calculate edit similarity
                    ident = 1.0 - (res["editDistance"] / len(t1))

                    if ident >= 0.75:
                        categories["method2_ontarget_75"] += 1
                    else:
                        categories["method2_nonspecific_75"] += 1

                    if ident >= 0.80:
                        categories["method2_ontarget_80"] += 1
                    else:
                        categories["method2_nonspecific_80"] += 1
                else:
                    # Missing reference amplicon CANNOT be verified as on-target
                    categories["method2_nonspecific_75"] += 1
                    categories["method2_nonspecific_80"] += 1
            else:
                categories["method1_nonspecific"] += 1
                categories["method2_nonspecific_75"] += 1
                categories["method2_nonspecific_80"] += 1

            if max_reads and total_reads >= max_reads:
                break

    elapsed = time.time() - t0
    total_classified = categories["dimers"] + categories["method1_ontarget"] + categories["method1_nonspecific"]
    total_excluded = categories["no_adapter"] + categories["qc_failed"] + categories["primer_unmatched"]

    print("-" * 78)
    print("FASTQ PROCESSING RESULTS:")
    print(f"Elapsed Time:               {elapsed:.2f} s ({total_reads / elapsed:.0f} reads/s)")
    print(f"Total Processed Reads:      {total_reads:,}")
    print(f"Total Classified Reads:     {total_classified:,} ({total_classified/total_reads*100:.2f}%)")
    print(f"Total Excluded Reads:       {total_excluded:,} ({total_excluded/total_reads*100:.2f}%)")
    print(f"  - Missing Adapter:        {categories['no_adapter']:,} ({categories['no_adapter']/total_reads*100:.2f}%)")
    print(f"  - Quality / Length Filter:{categories['qc_failed']:,} ({categories['qc_failed']/total_reads*100:.2f}%)")
    print(f"  - Primer Unmatched:       {categories['primer_unmatched']:,} ({categories['primer_unmatched']/total_reads*100:.2f}%)")

    dimer_pct = categories["dimers"] / total_classified * 100 if total_classified > 0 else 0
    m1_on_pct = categories["method1_ontarget"] / total_classified * 100 if total_classified > 0 else 0
    m1_non_pct = categories["method1_nonspecific"] / total_classified * 100 if total_classified > 0 else 0

    m2_75_on = categories["method2_ontarget_75"] / total_classified * 100 if total_classified > 0 else 0
    m2_75_non = categories["method2_nonspecific_75"] / total_classified * 100 if total_classified > 0 else 0

    m2_80_on = categories["method2_ontarget_80"] / total_classified * 100 if total_classified > 0 else 0
    m2_80_non = categories["method2_nonspecific_80"] / total_classified * 100 if total_classified > 0 else 0

    # Dynamic status calculation
    dimer_diff = abs(dimer_pct - 1.0)
    dimer_status = f"MATCH (diff = {dimer_pct - 1.0:+.2f}%)" if dimer_diff <= 0.5 else f"MISMATCH ({dimer_pct:.2f}% vs ~1%)"

    m1_diff = abs(m1_on_pct - 43.0)
    m1_status = f"MATCH (diff = {m1_on_pct - 43.0:+.2f}%)" if m1_diff <= 1.0 else f"MISMATCH ({m1_on_pct:.2f}% vs ~43%)"

    m2_75_diff = abs(m2_75_on - 43.0)
    m2_75_status = f"MATCH (diff = {m2_75_on - 43.0:+.2f}%)" if m2_75_diff <= 1.0 else f"MISMATCH ({m2_75_on:.2f}% vs ~43%)"

    m2_80_diff = abs(m2_80_on - 43.0)
    m2_80_status = f"MATCH (diff = {m2_80_on - 43.0:+.2f}%)" if m2_80_diff <= 1.0 else f"MISMATCH ({m2_80_on:.2f}% vs ~43%)"

    print("\n[Read Category Breakdown]:")
    print(f"  Primer-Dimers (invariant):     {categories['dimers']:,} ({dimer_pct:.2f}%)  | Published Fig 3f: ~1% [{dimer_status}]")
    print(f"\n  Method 1 (Primer Pairing):")
    print(f"    - On-Target Amplicons:       {categories['method1_ontarget']:,} ({m1_on_pct:.2f}%)  | Published Fig 3f: ~43% [{m1_status}]")
    print(f"    - Nonspecific Products:      {categories['method1_nonspecific']:,} ({m1_non_pct:.2f}%)  | Published Fig 3f: ~56%")
    print(f"\n  Method 2 (Cognate-Anchored Gapped Alignment via edlib >= 75% Identity):")
    print(f"    - On-Target Amplicons:       {categories['method2_ontarget_75']:,} ({m2_75_on:.2f}%)  | Published Fig 3f: ~43% [{m2_75_status}]")
    print(f"    - Nonspecific Products:      {categories['method2_nonspecific_75']:,} ({m2_75_non:.2f}%)  | Published Fig 3f: ~56%")
    print(f"\n  Method 2 (Cognate-Anchored Gapped Alignment via edlib >= 80% Identity):")
    print(f"    - On-Target Amplicons:       {categories['method2_ontarget_80']:,} ({m2_80_on:.2f}%)  | Published Fig 3f: ~43% [{m2_80_status}]")
    print(f"    - Nonspecific Products:      {categories['method2_nonspecific_80']:,} ({m2_80_non:.2f}%)  | Published Fig 3f: ~56%")

    # Dynamic pass/fail verification
    all_matched = (dimer_diff <= 0.5) and (m1_diff <= 1.0) and (m2_75_diff <= 1.0 or m2_80_diff <= 1.0)
    print("\n" + "=" * 78)
    if all_matched:
        print("STATUS: FASTQ CLASSIFICATION SUCCESSFULLY COMPLETED AND REPLICATED")
        print("=" * 78)
        return True
    else:
        print("STATUS: FASTQ CLASSIFICATION COMPLETED - REPLICATION TOLERANCE NOT MET", file=sys.stderr)
        print("=" * 78, file=sys.stderr)
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reproduce SADDLE Read Classification Directly from Raw FASTQ")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--full", action="store_true", help="Process all 3,942,385 reads in the E19 dataset")
    group.add_argument("--sample", type=int, nargs="?", const=50000, help="Process the first N reads (default: 50,000)")
    parser.add_argument("--raw-dir", type=str, default=None, help="Custom path to raw data directory")
    args = parser.parse_args()

    raw_directory = find_raw_dir(args.raw_dir)
    max_count = None if args.full else (args.sample if args.sample else 50000)
    success = run_pipeline(raw_directory, max_reads=max_count)
    sys.exit(0 if success else 1)
