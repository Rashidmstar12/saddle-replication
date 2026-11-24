#!/usr/bin/env python3
"""
code/run_saddle_replication.py

Open-source computational replication of SADDLE (Xie et al., Nature
Communications 13, 1881, 2022) across four primary NGS validation panels:
  - 96-plex PS1 (Fig. 3a; published primer-dimers: 90.7%)
  - 96-plex PS2 (Fig. 3a; published primer-dimers: 39.6%)
  - 96-plex PS3 (Fig. 3a; published primer-dimers: 4.9%)
  - 384-plex E19 (Fig. 3f; published: ~43% on-target, ~56% nonspecific, ~1% dimers)

V4 enhancements:
  - Fail-fast Bowtie2 executable check before reading raw FASTQ files
    (requires explicit --skip-bowtie2 for partial Pipeline B runs).
  - Archives all 6 reconstructed 96-plex matrices with both directional ordered
    counts (reads_p1_p2, reads_p2_p1) and unordered totals (dimer_reads),
    preserving the complete 192x192 ordered matrix and 18,528 unordered pairs.
  - Reports active-pair-union exact agreement and count-stratified disagreement
    breakdowns directly against published MOESM7 matrices.
"""

import argparse
import gzip
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import openpyxl
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import canonical_classifier as cc


PANELS = [
    {
        "panel_id": "PS1",
        "plex": 96,
        "figure_ref": "Figure 3a",
        "r1": "SADDLE_PS1_E14_Raw_R1.fastq.gz",
        "r2": "SADDLE_PS1_E14_Raw_R2.fastq.gz",
        "primer_sheet": "96plex PS1",
        "ref_fasta": "Ref_1.fasta",
        "moesm7_not_opt": "PS1_Dimer NGS(not optimized)",
        "moesm7_conc_opt": "PS1_Dimer NGS (conc optimized)",
        "published": {"dimer_pct": 90.7, "ontarget_pct": None, "nonspecific_pct": None},
    },
    {
        "panel_id": "PS2",
        "plex": 96,
        "figure_ref": "Figure 3a",
        "r1": "SADDLE_PS2_E14_Raw_R1.fastq.gz",
        "r2": "SADDLE_PS2_E14_Raw_R2.fastq.gz",
        "primer_sheet": "96plex PS2",
        "ref_fasta": "Ref_2.fasta",
        "moesm7_not_opt": "PS2_Dimer NGS(not optimized)",
        "moesm7_conc_opt": "PS2_Dimer NGS(conc optimized)",
        "published": {"dimer_pct": 39.6, "ontarget_pct": None, "nonspecific_pct": None},
    },
    {
        "panel_id": "PS3",
        "plex": 96,
        "figure_ref": "Figure 3a",
        "r1": "SADDLE_PS3_E14_Raw_R1.fastq.gz",
        "r2": "SADDLE_PS3_E14_Raw_R2.fastq.gz",
        "primer_sheet": "96plex PS3",
        "ref_fasta": "Ref_3.fasta",
        "moesm7_not_opt": "PS3_Dimer NGS(not optimized)",
        "moesm7_conc_opt": "PS3_Dimer NGS(conc optimized)",
        "published": {"dimer_pct": 4.9, "ontarget_pct": None, "nonspecific_pct": None},
    },
    {
        "panel_id": "E19",
        "plex": 384,
        "figure_ref": "Figure 3f",
        "r1": "E19_384plex_NGS_R1_fastq.gz",
        "r2": "E19_384plex_NGS_R2_fastq.gz",
        "primer_sheet": "384 Plex Primer",
        "ref_fasta": "reference_384_amplicons.fasta",
        "moesm7_not_opt": None,
        "moesm7_conc_opt": None,
        "published": {"dimer_pct": 1.0, "ontarget_pct": 43.0, "nonspecific_pct": 56.0},
    },
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_raw_input_hashes(
    input_files: List[Path], expected_manifest_path: Path, skip_check: bool = False
) -> dict:
    expected_map = {}
    if expected_manifest_path.exists():
        expected_map = json.loads(expected_manifest_path.read_text(encoding="utf-8")).get("files", {})

    provenance = {}
    for fp in input_files:
        if not fp.exists():
            raise FileNotFoundError(f"Missing required raw input file: {fp}")
        size = fp.stat().st_size
        digest = sha256_file(fp)
        verified = False
        if fp.name in expected_map:
            exp = expected_map[fp.name]
            if not skip_check:
                if size != exp["size_bytes"]:
                    raise ValueError(
                        f"Raw input size mismatch for {fp.name}: got {size}, expected {exp['size_bytes']}"
                    )
                if digest != exp["sha256"]:
                    raise ValueError(
                        f"Raw input SHA-256 mismatch for {fp.name}: got {digest}, expected {exp['sha256']}"
                    )
                verified = True
        elif not skip_check:
            raise ValueError(f"Raw input {fp.name} not registered in {expected_manifest_path.name}")
        provenance[fp.name] = {
            "size_bytes": size,
            "sha256": digest,
            "verified_against_expected_manifest": verified,
        }
    return provenance


def find_raw_dir(base_dir: Path, custom: Optional[str] = None) -> Path:
    if custom:
        p = Path(custom).resolve()
        if p.exists():
            return p
    for cand in [base_dir / "data" / "raw", base_dir.parent / "data" / "raw", Path("data/raw")]:
        if cand.exists() and (cand / "41467_2022_29500_MOESM5_ESM.xlsx").exists():
            return cand.resolve()
    raise FileNotFoundError("Could not locate raw data directory containing MOESM5 and FASTQs.")


def read_fasta_sequences(fasta_path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name = None
    chunks: List[str] = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).lower()))
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).lower()))
    return records


def extract_strict_amplicon(fp_name: str, fp_seq: str, rp_seq: str, full_seq: str) -> str:
    fp_low = fp_seq.strip().lower()
    rp_rc = cc.revcomp(rp_seq.strip().lower())
    ctx = full_seq.strip().lower()
    s_pos = ctx.find(fp_low)
    if s_pos == -1:
        raise ValueError(f"Forward primer {fp_name} not found in genomic context sequence.")
    e_pos = ctx.find(rp_rc, s_pos)
    if e_pos == -1:
        raise ValueError(f"Reverse-complement reverse primer for {fp_name} not found after forward primer.")
    amp = ctx[s_pos : e_pos + len(rp_rc)]
    if not (amp.startswith(fp_low) and amp.endswith(rp_rc)):
        raise ValueError(f"Extracted amplicon for {fp_name} failed boundary assertion.")
    return amp


def build_384_reference_fasta(moesm5_path: Path, out_fasta: Path) -> List[str]:
    wb = openpyxl.load_workbook(moesm5_path, read_only=True, data_only=True)
    ws = wb["384 Plex Info"]
    rows = list(ws.iter_rows(values_only=True))[1:]
    wb.close()
    if len(rows) != 384:
        raise ValueError(f"Expected 384 amplicon rows in MOESM5 '384 Plex Info', got {len(rows)}")

    amplicons: List[str] = []
    lines: List[str] = []
    for idx, r in enumerate(rows, 1):
        fp_name = str(r[0]).strip()
        amp = extract_strict_amplicon(fp_name, str(r[1]), str(r[3]), str(r[9]))
        amplicons.append(amp)
        lines.append(f">amplicon_{idx}_{fp_name}\n{amp.upper()}\n")

    cc.write_utf8_lf(out_fasta, "".join(lines))
    return amplicons


def load_panel_primers_and_refs(
    panel_cfg: dict, moesm5_path: Path, source_ref_dir: Path, out_ref_dir: Path
) -> Tuple[List[Tuple[str, str]], List[str], Dict[str, List[int]], List[str], Path]:
    wb = openpyxl.load_workbook(moesm5_path, read_only=True, data_only=True)
    ws = wb[panel_cfg["primer_sheet"]]
    fps: List[Tuple[str, str]] = []
    rps: List[Tuple[str, str]] = []
    for r in ws.iter_rows(values_only=True):
        if r[0] is None:
            continue
        fps.append((str(r[0]).strip(), str(r[1]).strip().lower()))
        rps.append((str(r[2]).strip(), str(r[3]).strip().lower()))
    wb.close()

    plex = panel_cfg["plex"]
    if len(fps) != plex or len(rps) != plex:
        raise ValueError(f"{panel_cfg['panel_id']}: expected {plex} primer pairs, got ({len(fps)}, {len(rps)})")
    all_primers = fps + rps
    primers_lower = [s for _, s in all_primers]
    prefix_index = cc.build_prefix_index(all_primers, k=8)

    if panel_cfg["panel_id"] == "E19":
        ref_fasta = out_ref_dir / panel_cfg["ref_fasta"]
        amplicons = build_384_reference_fasta(moesm5_path, ref_fasta)
    else:
        ref_fasta = source_ref_dir / panel_cfg["ref_fasta"]
        recs = read_fasta_sequences(ref_fasta)
        if len(recs) != plex:
            raise ValueError(f"{ref_fasta}: expected {plex} seqs, got {len(recs)}")
        amplicons = [s for _, s in recs]
        for i in range(plex):
            fp_seq = fps[i][1]
            rp_rc = cc.revcomp(rps[i][1])
            if not (amplicons[i].startswith(fp_seq) and amplicons[i].endswith(rp_rc)):
                raise ValueError(f"{panel_cfg['panel_id']} amplicon {i + 1} failed primer boundary verification")

    return all_primers, primers_lower, prefix_index, amplicons, ref_fasta


def get_bowtie2_paths() -> Tuple[Optional[str], Optional[str], Optional[str]]:
    bt2 = shutil.which("bowtie2")
    bt2_build = shutil.which("bowtie2-build")
    if not bt2 or not bt2_build:
        return None, None, None
    try:
        res = subprocess.run([bt2, "--version"], capture_output=True, text=True, check=True)
        m = re.search(r"version\s+([0-9]+\.[0-9]+\.[0-9]+[^\s]*)", res.stdout or "")
        ver_line = f"bowtie2 {m.group(1)}" if m else "bowtie2 unknown"
    except Exception:
        ver_line = "bowtie2 unknown"
    return bt2, bt2_build, ver_line


def require_bowtie2(skip_bowtie2: bool) -> Tuple[bool, Optional[str], Optional[str], Optional[str]]:
    """Fail-fast dependency check: require Bowtie2 unless --skip-bowtie2 is explicitly passed."""
    if skip_bowtie2:
        return False, None, None, None
    bt2, bt2_build, bt2_ver = get_bowtie2_paths()
    if bt2 is None or bt2_build is None:
        raise RuntimeError(
            "Required executables 'bowtie2' and 'bowtie2-build' were not found on PATH. "
            "Install Bowtie2 (see environment.yml / README.md) or pass --skip-bowtie2 "
            "for a partial Pipeline B run."
        )
    return True, bt2, bt2_build, bt2_ver


def summarize_dimer_matrix(
    mat_ordered: np.ndarray, all_primers: List[Tuple[str, str]]
) -> dict:
    n = mat_ordered.shape[0]
    triu_i, triu_j = np.triu_indices(n)
    unordered = mat_ordered[triu_i, triu_j] + np.where(
        triu_i != triu_j, mat_ordered[triu_j, triu_i], 0
    )
    total_pairs = int(len(unordered))
    total_dimer_reads = int(unordered.sum())
    active_pairs = int((unordered > 0).sum())
    severe_pairs = int((unordered >= 10).sum())
    zero_pairs = total_pairs - active_pairs

    order = np.argsort(-unordered, kind="stable")
    top5 = []
    for rank_idx in range(min(5, total_pairs)):
        k = int(order[rank_idx])
        cnt = int(unordered[k])
        i, j = int(triu_i[k]), int(triu_j[k])
        top5.append({
            "rank": rank_idx + 1,
            "pair_index": k,
            "primer_1": all_primers[i][0],
            "primer_2": all_primers[j][0],
            "dimer_reads": cnt,
            "pct_of_dimer_reads": round((cnt / total_dimer_reads) * 100.0, 4) if total_dimer_reads else 0.0,
        })

    top10_reads = int(unordered[order[:10]].sum()) if total_dimer_reads else 0
    return {
        "total_candidate_pairs": total_pairs,
        "total_dimer_reads": total_dimer_reads,
        "zero_dimer_pairs": zero_pairs,
        "zero_dimer_pct": round((zero_pairs / total_pairs) * 100.0, 4),
        "active_dimer_pairs_ge1": active_pairs,
        "active_dimer_pct_ge1": round((active_pairs / total_pairs) * 100.0, 4),
        "severe_dimer_pairs_ge10": severe_pairs,
        "severe_dimer_pct_ge10": round((severe_pairs / total_pairs) * 100.0, 4),
        "top1_pct_of_dimer_reads": top5[0]["pct_of_dimer_reads"] if top5 else 0.0,
        "top10_pct_of_dimer_reads": round((top10_reads / total_dimer_reads) * 100.0, 4) if total_dimer_reads else 0.0,
        "top5_pairs": top5,
    }


def write_pair_table_tsv(
    mat_ordered: np.ndarray,
    all_primers: List[Tuple[str, str]],
    out_tsv: Path,
    include_directional: bool = False,
    compress_gz: bool = False,
) -> None:
    """
    Write pair table with deterministic LF endings (and mtime=0 for .gz).
    When include_directional=True, writes reads_p1_p2 and reads_p2_p1 alongside
    dimer_reads so the full ordered 2N x 2N matrix and unordered pair totals are
    both preserved losslessly in one compact file.
    """
    n = mat_ordered.shape[0]
    triu_i, triu_j = np.triu_indices(n)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    if include_directional:
        lines = ["pair_index\tprimer_1_id\tprimer_2_id\treads_p1_p2\treads_p2_p1\tdimer_reads\n"]
        for k in range(len(triu_i)):
            i, j = int(triu_i[k]), int(triu_j[k])
            fwd = int(mat_ordered[i, j])
            rev = int(mat_ordered[j, i]) if i != j else 0
            lines.append(f"{k}\t{all_primers[i][0]}\t{all_primers[j][0]}\t{fwd}\t{rev}\t{fwd + rev}\n")
    else:
        lines = ["pair_index\tprimer_1_id\tprimer_2_id\tdimer_reads\n"]
        for k in range(len(triu_i)):
            i, j = int(triu_i[k]), int(triu_j[k])
            tot = int(mat_ordered[i, j]) + (int(mat_ordered[j, i]) if i != j else 0)
            lines.append(f"{k}\t{all_primers[i][0]}\t{all_primers[j][0]}\t{tot}\n")

    payload = "".join(lines).encode("utf-8")
    if compress_gz:
        with open(out_tsv, "wb") as raw_f:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_f, mtime=0) as gz:
                gz.write(payload)
    else:
        with open(out_tsv, "wb") as raw_f:
            raw_f.write(payload)


def compare_with_moesm7(
    mat_ordered: np.ndarray,
    all_primers: List[Tuple[str, str]],
    moesm7_path: Path,
    sheet_not_opt: str,
    sheet_conc_opt: str,
) -> dict:
    """
    Compare reconstructed 192x192 dimer matrix against Xie et al. MOESM7 with:
      - Pearson r (unordered & ordered), Spearman rho
      - Exact-cell agreement (all 18,528 pairs AND restricted to active-pair union)
      - Active-pair Jaccard overlap, MAE (all pairs & active union), normalized L1 error
      - Count-stratified disagreement breakdown (max(rec, pub) in 1-2, 3-9, >=10 reads)
    """
    wb7 = openpyxl.load_workbook(moesm7_path, read_only=True, data_only=True)
    n = mat_ordered.shape[0]
    triu_i, triu_j = np.triu_indices(n)
    rec_unord = (
        mat_ordered[triu_i, triu_j] + np.where(triu_i != triu_j, mat_ordered[triu_j, triu_i], 0)
    ).astype(np.int64)
    rec_ord = mat_ordered.ravel().astype(np.int64)

    out = {}
    for label, sname in [("not_optimized", sheet_not_opt), ("conc_optimized", sheet_conc_opt)]:
        pub_mat = np.array(
            [[int(round(float(c or 0))) for c in r] for r in wb7[sname].iter_rows(values_only=True)],
            dtype=np.int64,
        )
        pub_unord = (
            pub_mat[triu_i, triu_j] + np.where(triu_i != triu_j, pub_mat[triu_j, triu_i], 0)
        ).astype(np.int64)
        pub_ord = pub_mat.ravel().astype(np.int64)

        r_unord = float(pearsonr(rec_unord, pub_unord)[0]) if rec_unord.std() > 0 and pub_unord.std() > 0 else 0.0
        rho_unord = float(spearmanr(rec_unord, pub_unord)[0])
        r_ord = float(pearsonr(rec_ord, pub_ord)[0]) if rec_ord.std() > 0 and pub_ord.std() > 0 else 0.0

        abs_diff = np.abs(rec_unord - pub_unord)
        exact_all_pct = float(np.mean(rec_unord == pub_unord) * 100.0)
        rec_active = rec_unord > 0
        pub_active = pub_unord > 0
        union_active = rec_active | pub_active
        inter_active = rec_active & pub_active
        n_union = int(union_active.sum())
        exact_active_union_count = int(np.sum((rec_unord == pub_unord) & union_active))
        exact_active_union_pct = float((exact_active_union_count / n_union) * 100.0) if n_union > 0 else 100.0
        jaccard = float(inter_active.sum() / n_union) if n_union > 0 else 1.0
        mae_all = float(abs_diff.mean())
        mae_union = float(abs_diff[union_active].mean()) if n_union > 0 else 0.0
        total_l1 = int(abs_diff.sum())
        norm_l1 = float(total_l1 / pub_unord.sum()) if pub_unord.sum() > 0 else 0.0

        # Count-stratified disagreement breakdown (directly verifiable!)
        disagree = rec_unord != pub_unord
        n_dis = int(disagree.sum())
        max_cnt = np.maximum(rec_unord, pub_unord)
        low_mask = disagree & (max_cnt <= 2)
        mid_mask = disagree & (max_cnt >= 3) & (max_cnt <= 9)
        high_mask = disagree & (max_cnt >= 10)

        def strat_entry(mask: np.ndarray) -> dict:
            cnt = int(mask.sum())
            l1 = int(abs_diff[mask].sum())
            return {
                "cell_count": cnt,
                "pct_of_disagreeing_cells": round((cnt / n_dis) * 100.0, 2) if n_dis else 0.0,
                "l1_reads": l1,
                "pct_of_total_l1": round((l1 / total_l1) * 100.0, 2) if total_l1 else 0.0,
            }

        rec_top10 = set(int(x) for x in np.argsort(-rec_unord, kind="stable")[:10])
        pub_top10 = set(int(x) for x in np.argsort(-pub_unord, kind="stable")[:10])
        rec_top1 = int(np.argsort(-rec_unord, kind="stable")[0])
        pub_top1 = int(np.argsort(-pub_unord, kind="stable")[0])

        out[label] = {
            "sheet_name": sname,
            "published_total_dimer_reads": int(pub_mat.sum()),
            "reconstructed_total_dimer_reads": int(mat_ordered.sum()),
            "total_read_ratio_rec_over_pub": round(float(mat_ordered.sum()) / float(pub_mat.sum()), 4)
            if pub_mat.sum() > 0
            else None,
            "pearson_r_unordered": round(r_unord, 6),
            "spearman_rho_unordered": round(rho_unord, 6),
            "pearson_r_ordered": round(r_ord, 6),
            "exact_cell_match_pct": round(exact_all_pct, 4),
            "exact_cell_match_active_union_count": exact_active_union_count,
            "exact_cell_match_active_union_pct": round(exact_active_union_pct, 2),
            "active_pairs_union": n_union,
            "active_pair_jaccard": round(jaccard, 4),
            "active_pairs_reconstructed": int(rec_active.sum()),
            "active_pairs_published": int(pub_active.sum()),
            "active_pairs_intersection": int(inter_active.sum()),
            "mae_all_pairs": round(mae_all, 4),
            "mae_active_union": round(mae_union, 4),
            "total_l1_error_reads": total_l1,
            "normalized_l1_error": round(norm_l1, 6),
            "disagreeing_cells_total": n_dis,
            "disagreement_stratification": {
                "max_count_1_to_2": strat_entry(low_mask),
                "max_count_3_to_9": strat_entry(mid_mask),
                "max_count_ge_10": strat_entry(high_mask),
            },
            "top1_exact_match": bool(rec_top1 == pub_top1),
            "top10_overlap_count": int(len(rec_top10 & pub_top10)),
            "published_top1_pair": f"{all_primers[triu_i[pub_top1]][0]} x {all_primers[triu_j[pub_top1]][0]}",
            "reconstructed_top1_pair": f"{all_primers[triu_i[rec_top1]][0]} x {all_primers[triu_j[rec_top1]][0]}",
        }
    wb7.close()
    return out


def run_panel_pipeline_a_bowtie2(
    panel_cfg: dict,
    r1_path: Path,
    r2_path: Path,
    all_primers: List[Tuple[str, str]],
    primers_lower: List[str],
    prefix_index: Dict[str, List[int]],
    ref_fasta: Path,
    bt2: str,
    bt2_build: str,
    threads: int = 4,
    max_reads: Optional[int] = None,
) -> Tuple[dict, np.ndarray]:
    t0 = time.time()
    n_primers = len(primers_lower)
    dimer_mat = np.zeros((n_primers, n_primers), dtype=np.int64)
    others_mat = np.zeros((n_primers, n_primers), dtype=np.int64)

    counts = {
        "total_raw_reads": 0,
        "no_adapter_both": 0,
        "no_adapter_r1_only": 0,
        "no_adapter_r2_only": 0,
        "qc_short_len": 0,
        "qc_ambiguous_n": 0,
        "qc_low_quality": 0,
        "qc_passed": 0,
        "routed_long_ge60_to_bowtie2": 0,
        "routed_short_lt60_to_dimer_cand": 0,
        "mixed_length_60bp_dropped_matlab_rule": 0,
        "bt2_both_mapped_concordant": 0,
        "bt2_both_mapped_discordant": 0,
        "bt2_one_mate_mapped": 0,
        "bt2_both_unmapped": 0,
        "bt2_mapped_reads_total": 0,
        "dimer_cand_total": 0,
        "dimer_cand_len_le15_excluded": 0,
        "dimer_cand_primer_unmatched": 0,
        "dimers": 0,
        "nonspecific": 0,
    }

    with tempfile.TemporaryDirectory(prefix=f"saddle_bt2_{panel_cfg['panel_id']}_") as tmp_str:
        tmp = Path(tmp_str)
        long_r1 = tmp / "long_R1.fastq"
        long_r2 = tmp / "long_R2.fastq"
        short_pairs: List[Tuple[str, str]] = []

        with (
            open(long_r1, "w", encoding="ascii", newline="\n") as out1,
            open(long_r2, "w", encoding="ascii", newline="\n") as out2,
        ):
            for _, s1, q1, s2, q2 in cc.iter_validated_fastq_pairs(r1_path, r2_path, max_reads=max_reads):
                counts["total_raw_reads"] += 1
                t1, tq1, t2, tq2, adp_st = cc.trim_adapter_pair(
                    s1, q1, s2, q2, cc.ORIG_ADAPTER_R1, cc.ORIG_ADAPTER_R2
                )
                if adp_st != "ok":
                    counts[adp_st] += 1
                    continue

                qc_st = cc.qc_check_pair(t1, tq1, t2, tq2, min_len_exclusive=10, q30_frac_exclusive=0.7)
                if qc_st != "ok":
                    counts[qc_st] += 1
                    continue

                counts["qc_passed"] += 1
                l1, l2 = len(t1), len(t2)
                if l1 >= 60 and l2 >= 60:
                    counts["routed_long_ge60_to_bowtie2"] += 1
                    rname = f"@r{counts['total_raw_reads']}"
                    out1.write(f"{rname}\n{t1}\n+\n{tq1}\n")
                    out2.write(f"{rname}\n{t2}\n+\n{tq2}\n")
                elif l1 < 60 and l2 < 60:
                    counts["routed_short_lt60_to_dimer_cand"] += 1
                    short_pairs.append((t1.lower(), t2.lower()))
                else:
                    counts["mixed_length_60bp_dropped_matlab_rule"] += 1

        unmapped_long_pairs: List[Tuple[str, str]] = []
        if counts["routed_long_ge60_to_bowtie2"] > 0:
            idx_prefix = tmp / "ref_idx"
            subprocess.run(
                [bt2_build, "--quiet", str(ref_fasta), str(idx_prefix)],
                check=True,
                capture_output=True,
                text=True,
            )
            sam_path = tmp / "aligned.sam"
            bt2_proc = subprocess.run(
                [
                    bt2,
                    "-x",
                    str(idx_prefix),
                    "-1",
                    str(long_r1),
                    "-2",
                    str(long_r2),
                    "--threads",
                    str(threads),
                    "-S",
                    str(sam_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            bt2_stderr = bt2_proc.stderr.strip()

            pending: Dict[str, Tuple[int, str]] = {}
            with open(sam_path, "r", encoding="ascii") as sf:
                for line in sf:
                    if not line or line[0] == "@":
                        continue
                    parts = line.split("\t", 11)
                    qname = parts[0]
                    flag = int(parts[1])
                    if flag & 0x900:
                        continue
                    seq = parts[9]
                    if qname not in pending:
                        pending[qname] = (flag, seq)
                    else:
                        f_prev, s_prev = pending.pop(qname)
                        if f_prev & 0x40:
                            f_r1, s_r1, f_r2, s_r2 = f_prev, s_prev, flag, seq
                        else:
                            f_r1, s_r1, f_r2, s_r2 = flag, seq, f_prev, s_prev

                        r1_unm = bool(f_r1 & 0x4)
                        r2_unm = bool(f_r2 & 0x4)
                        if not r1_unm:
                            counts["bt2_mapped_reads_total"] += 1
                        if not r2_unm:
                            counts["bt2_mapped_reads_total"] += 1

                        if (not r1_unm) and (not r2_unm):
                            if f_r1 & 0x2:
                                counts["bt2_both_mapped_concordant"] += 1
                            else:
                                counts["bt2_both_mapped_discordant"] += 1
                        elif r1_unm and r2_unm:
                            counts["bt2_both_unmapped"] += 1
                            unmapped_long_pairs.append((s_r1.lower(), s_r2.lower()))
                        else:
                            counts["bt2_one_mate_mapped"] += 1
        else:
            bt2_stderr = ""

        for t1_low, t2_low in short_pairs + unmapped_long_pairs:
            counts["dimer_cand_total"] += 1
            if len(t1_low) <= 15:
                counts["dimer_cand_len_le15_excluded"] += 1
                continue
            cands1 = prefix_index.get(t1_low[:8])
            cands2 = prefix_index.get(t2_low[:8])
            if not cands1 or not cands2:
                counts["dimer_cand_primer_unmatched"] += 1
                continue
            id1 = cc.match_primer_sw_nuc44(t1_low, cands1, primers_lower, threshold=0.6)
            if id1 is None:
                counts["dimer_cand_primer_unmatched"] += 1
                continue
            id2 = cc.match_primer_sw_nuc44(t2_low, cands2, primers_lower, threshold=0.6)
            if id2 is None:
                counts["dimer_cand_primer_unmatched"] += 1
                continue

            if cc.is_dimer_insert(len(t1_low), len(primers_lower[id1]), len(primers_lower[id2])):
                counts["dimers"] += 1
                dimer_mat[id1, id2] += 1
            else:
                counts["nonspecific"] += 1
                others_mat[id1, id2] += 1

    ontarget_ipynb_formula = counts["bt2_mapped_reads_total"] // 2
    ontarget_both_mapped = counts["bt2_both_mapped_concordant"] + counts["bt2_both_mapped_discordant"]

    classified_ipynb = ontarget_ipynb_formula + counts["dimers"] + counts["nonspecific"]
    classified_both = ontarget_both_mapped + counts["dimers"] + counts["nonspecific"]
    raw_total = counts["total_raw_reads"]

    def pct(num: int, den: int) -> float:
        return round((num / den) * 100.0, 4) if den > 0 else 0.0

    elapsed = round(time.time() - t0, 2)
    result = {
        "elapsed_seconds": elapsed,
        "counts": counts,
        "ontarget_ipynb_formula_pairs": ontarget_ipynb_formula,
        "ontarget_both_mapped_pairs": ontarget_both_mapped,
        "ontarget_floor_div_extra_pairs": ontarget_ipynb_formula - ontarget_both_mapped,
        "classified_denominator_ipynb": classified_ipynb,
        "classified_denominator_both_mapped": classified_both,
        "pct_of_classified_ipynb": {
            "ontarget_pct": pct(ontarget_ipynb_formula, classified_ipynb),
            "nonspecific_pct": pct(counts["nonspecific"], classified_ipynb),
            "dimer_pct": pct(counts["dimers"], classified_ipynb),
        },
        "pct_of_classified_both_mapped": {
            "ontarget_pct": pct(ontarget_both_mapped, classified_both),
            "nonspecific_pct": pct(counts["nonspecific"], classified_both),
            "dimer_pct": pct(counts["dimers"], classified_both),
        },
        "pct_of_raw_reads": {
            "ontarget_ipynb_pct": pct(ontarget_ipynb_formula, raw_total),
            "nonspecific_pct": pct(counts["nonspecific"], raw_total),
            "dimer_pct": pct(counts["dimers"], raw_total),
            "classified_ipynb_pct": pct(classified_ipynb, raw_total),
            "excluded_pct": pct(raw_total - classified_ipynb, raw_total),
        },
        "dimer_distribution": summarize_dimer_matrix(dimer_mat, all_primers),
        "bowtie2_summary_stderr": bt2_stderr,
    }
    return result, dimer_mat


def run_panel_pipeline_b_edlib(
    panel_cfg: dict,
    r1_path: Path,
    r2_path: Path,
    all_primers: List[Tuple[str, str]],
    primers_lower: List[str],
    prefix_index: Dict[str, List[int]],
    amplicons: List[str],
    max_reads: Optional[int] = None,
) -> Tuple[dict, np.ndarray]:
    t0 = time.time()
    plex = panel_cfg["plex"]
    n_primers = len(primers_lower)
    dimer_mat = np.zeros((n_primers, n_primers), dtype=np.int64)

    counts = {
        "total_raw_reads": 0,
        "no_adapter_both": 0,
        "no_adapter_r1_only": 0,
        "no_adapter_r2_only": 0,
        "qc_short_len": 0,
        "qc_ambiguous_n": 0,
        "qc_low_quality": 0,
        "primer_unmatched_both": 0,
        "primer_unmatched_r1_only": 0,
        "primer_unmatched_r2_only": 0,
        "classified_total": 0,
        "dimers": 0,
        "m1_ontarget": 0,
        "m1_nonspecific": 0,
        "m2_edlib75_ontarget": 0,
        "m2_edlib80_ontarget": 0,
        "m2_edlib85_ontarget": 0,
        "ungapped75_ontarget": 0,
    }

    for _, s1, q1, s2, q2 in cc.iter_validated_fastq_pairs(r1_path, r2_path, max_reads=max_reads):
        counts["total_raw_reads"] += 1
        t1, tq1, t2, tq2, adp_st = cc.trim_adapter_pair(
            s1, q1, s2, q2, cc.CORE_ADAPTER_R1, cc.CORE_ADAPTER_R2
        )
        if adp_st != "ok":
            counts[adp_st] += 1
            continue

        qc_st = cc.qc_check_pair(t1, tq1, t2, tq2, min_len_exclusive=10, q30_frac_exclusive=0.7)
        if qc_st != "ok":
            counts[qc_st] += 1
            continue

        t1_low = t1.lower()
        t2_low = t2.lower()
        id1 = cc.match_primer_prefix(t1_low, prefix_index.get(t1_low[:8], []), primers_lower)
        id2 = cc.match_primer_prefix(t2_low, prefix_index.get(t2_low[:8], []), primers_lower)

        if id1 is None and id2 is None:
            counts["primer_unmatched_both"] += 1
            continue
        if id1 is None:
            counts["primer_unmatched_r1_only"] += 1
            continue
        if id2 is None:
            counts["primer_unmatched_r2_only"] += 1
            continue

        counts["classified_total"] += 1
        if cc.is_dimer_insert(len(t1_low), len(primers_lower[id1]), len(primers_lower[id2])):
            counts["dimers"] += 1
            dimer_mat[id1, id2] += 1
        else:
            amp_idx: Optional[int] = None
            orient_read: Optional[str] = None
            if id1 < plex and id2 == id1 + plex:
                amp_idx = id1
                orient_read = t1_low
            elif id2 < plex and id1 == id2 + plex:
                amp_idx = id2
                orient_read = t2_low

            if amp_idx is not None:
                counts["m1_ontarget"] += 1
                ref_seq = amplicons[amp_idx]
                ident = cc.edlib_shw_identity(orient_read, ref_seq)
                if ident >= 0.75:
                    counts["m2_edlib75_ontarget"] += 1
                    if ident >= 0.80:
                        counts["m2_edlib80_ontarget"] += 1
                        if ident >= 0.85:
                            counts["m2_edlib85_ontarget"] += 1
                if cc.ungapped_positional_identity(orient_read, ref_seq) >= 0.75:
                    counts["ungapped75_ontarget"] += 1
            else:
                counts["m1_nonspecific"] += 1

    cls = counts["classified_total"]
    raw_total = counts["total_raw_reads"]
    dim = counts["dimers"]
    missing_adp = counts["no_adapter_both"] + counts["no_adapter_r1_only"] + counts["no_adapter_r2_only"]
    qc_excl = counts["qc_short_len"] + counts["qc_ambiguous_n"] + counts["qc_low_quality"]
    unm_excl = counts["primer_unmatched_both"] + counts["primer_unmatched_r1_only"] + counts["primer_unmatched_r2_only"]
    excl_total = missing_adp + qc_excl + unm_excl
    assert cls + excl_total == raw_total, "Accounting closure failed!"

    def pct(num: int, den: int) -> float:
        return round((num / den) * 100.0, 4) if den > 0 else 0.0

    def method_block(ontarget_n: int) -> dict:
        nonspec_n = cls - dim - ontarget_n
        return {
            "ontarget_count": ontarget_n,
            "ontarget_pct_classified": pct(ontarget_n, cls),
            "ontarget_pct_raw": pct(ontarget_n, raw_total),
            "nonspecific_count": nonspec_n,
            "nonspecific_pct_classified": pct(nonspec_n, cls),
            "nonspecific_pct_raw": pct(nonspec_n, raw_total),
            "dimer_count": dim,
            "dimer_pct_classified": pct(dim, cls),
            "dimer_pct_raw": pct(dim, raw_total),
        }

    elapsed = round(time.time() - t0, 2)
    result = {
        "elapsed_seconds": elapsed,
        "counts": counts,
        "exclusions_summary": {
            "total_excluded": excl_total,
            "total_excluded_pct_raw": pct(excl_total, raw_total),
            "missing_adapter": missing_adp,
            "missing_adapter_pct_raw": pct(missing_adp, raw_total),
            "qc_failed": qc_excl,
            "qc_failed_pct_raw": pct(qc_excl, raw_total),
            "primer_unmatched": unm_excl,
            "primer_unmatched_pct_raw": pct(unm_excl, raw_total),
            "classified_pct_raw": pct(cls, raw_total),
        },
        "method1_primer_pairing": method_block(counts["m1_ontarget"]),
        "method2_edlib_75": method_block(counts["m2_edlib75_ontarget"]),
        "method2_edlib_80": method_block(counts["m2_edlib80_ontarget"]),
        "method2_edlib_85": method_block(counts["m2_edlib85_ontarget"]),
        "ungapped_positional_75": method_block(counts["ungapped75_ontarget"]),
        "dimer_distribution": summarize_dimer_matrix(dimer_mat, all_primers),
    }
    return result, dimer_mat


def resolve_ref_dir(base_dir: Path, raw_dir: Path, cli_ref_dir: Optional[str]) -> Path:
    candidates: List[Path] = []
    if cli_ref_dir:
        p = Path(cli_ref_dir).resolve()
        if (p / "Ref_1.fasta").exists():
            return p
        if (p / "reference_fastas" / "Ref_1.fasta").exists():
            return p / "reference_fastas"
        raise FileNotFoundError(f"Ref_1.fasta not found under --ref-dir {p}")
    candidates.extend([
        raw_dir / "reference_fastas",
        raw_dir,
        base_dir / "data" / "reference_fastas",
    ])
    for c in candidates:
        if (c / "Ref_1.fasta").exists():
            return c
    raise FileNotFoundError(
        "96-plex reference FASTAs (Ref_1.fasta-Ref_3.fasta) were not found. "
        "Fetch them programmatically from the pinned upstream commit via:\n"
        f"  python code/download_raw_data.py --refs-only --out-dir \"{raw_dir}\""
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run full multi-panel SADDLE open-source replication.")
    parser.add_argument("--raw-dir", type=str, default=None, help="Directory with raw FASTQ and MOESM files")
    parser.add_argument("--ref-dir", type=str, default=None, help="Directory with pinned Ref_1-3.fasta files (default: <raw-dir>/reference_fastas)")
    parser.add_argument("--out-dir", type=str, default=None, help="Output data directory (default: package data/)")
    parser.add_argument("--panels", nargs="+", default=["PS1", "PS2", "PS3", "E19"], help="Panels to execute")
    parser.add_argument("--threads", type=int, default=4, help="Bowtie2 threads")
    parser.add_argument("--max-reads", type=int, default=None, help="Optional cap for smoke testing")
    parser.add_argument("--skip-bowtie2", action="store_true", help="Explicitly skip Pipeline A (partial Pipeline B run)")
    parser.add_argument("--skip-hash-check", action="store_true", help="Skip raw SHA-256 verification")
    args = parser.parse_args()

    # Fail-fast Bowtie2 check BEFORE touching raw files
    run_bt2, bt2, bt2_build, bt2_ver = require_bowtie2(args.skip_bowtie2)

    base_dir = Path(__file__).resolve().parent.parent
    source_data_dir = base_dir / "data"
    out_data_dir = Path(args.out_dir).resolve() if args.out_dir else source_data_dir
    out_data_dir.mkdir(parents=True, exist_ok=True)
    out_ref_dir = out_data_dir / "reference_fastas"
    out_ref_dir.mkdir(parents=True, exist_ok=True)
    mat_dir = out_data_dir / "matrices"
    mat_dir.mkdir(parents=True, exist_ok=True)

    raw_dir = find_raw_dir(base_dir, args.raw_dir)
    selected = [p for p in PANELS if p["panel_id"] in args.panels]
    needs_96_refs = any(p["panel_id"] != "E19" for p in selected)
    source_ref_dir = resolve_ref_dir(base_dir, raw_dir, args.ref_dir) if needs_96_refs else (source_data_dir / "reference_fastas")
    moesm5_path = raw_dir / "41467_2022_29500_MOESM5_ESM.xlsx"
    moesm7_path = raw_dir / "41467_2022_29500_MOESM7_ESM.xlsx"

    print("=" * 78, flush=True)
    print("SADDLE MULTI-PANEL OPEN-SOURCE REPLICATION (PS1, PS2, PS3, E19)", flush=True)
    print(f"Raw dir:  {raw_dir}", flush=True)
    print(f"Ref dir:  {source_ref_dir}", flush=True)
    print(f"Out dir:  {out_data_dir}", flush=True)
    print(f"Bowtie2:  {bt2_ver if run_bt2 else 'EXPLICITLY SKIPPED (--skip-bowtie2)'}", flush=True)
    print("=" * 78, flush=True)

    input_files = [moesm5_path, moesm7_path]
    for p in selected:
        input_files.append(raw_dir / p["r1"])
        input_files.append(raw_dir / p["r2"])
        if p["panel_id"] != "E19":
            input_files.append(source_ref_dir / p["ref_fasta"])

    print("Verifying immutable raw input & pinned reference SHA-256 checksums...", flush=True)
    provenance = verify_raw_input_hashes(
        input_files,
        source_data_dir / "expected_raw_sha256.json",
        skip_check=args.skip_hash_check,
    )
    print(f"  Verified {len(provenance)} external input files.", flush=True)

    canonical: dict = {
        "schema_version": 5,
        "status": "DRAFT_UNAPPROVED",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "openpyxl": openpyxl.__version__,
            "bowtie2": bt2_ver if run_bt2 else None,
        },
        "raw_input_provenance": provenance,
        "panels": {},
    }

    for p in selected:
        pid = p["panel_id"]
        print(f"\n>>> Panel {pid} ({p['plex']}-plex, {p['figure_ref']}) <<<", flush=True)
        all_primers, primers_lower, prefix_index, amplicons, ref_fasta = load_panel_primers_and_refs(
            p, moesm5_path, source_ref_dir, out_ref_dir
        )
        r1_path = raw_dir / p["r1"]
        r2_path = raw_dir / p["r2"]

        panel_entry: dict = {
            "panel_id": pid,
            "plex": p["plex"],
            "n_primers": len(all_primers),
            "figure_ref": p["figure_ref"],
            "reference_fasta": ref_fasta.name,
            "reference_fasta_sha256": sha256_file(ref_fasta),
            "published_reference": p["published"],
        }

        if run_bt2:
            print(f"  [Pipeline A] Running Python + Bowtie2 on {pid}...", flush=True)
            res_a, mat_a = run_panel_pipeline_a_bowtie2(
                p, r1_path, r2_path, all_primers, primers_lower, prefix_index,
                ref_fasta, bt2, bt2_build, threads=args.threads, max_reads=args.max_reads,
            )
            if p["moesm7_not_opt"] is not None:
                res_a["moesm7_matrix_validation"] = compare_with_moesm7(
                    mat_a, all_primers, moesm7_path, p["moesm7_not_opt"], p["moesm7_conc_opt"]
                )
                if args.max_reads is None:
                    write_pair_table_tsv(
                        mat_a, all_primers, mat_dir / f"{pid.lower()}_pipeline_a_dimer_pairs.tsv.gz",
                        include_directional=True, compress_gz=True,
                    )
            elif pid == "E19" and args.max_reads is None:
                write_pair_table_tsv(
                    mat_a, all_primers, out_data_dir / "saddle_384_dimer_pairs_bowtie2.tsv",
                    include_directional=False, compress_gz=False,
                )

            panel_entry["pipeline_a_bowtie2"] = res_a
            pa_pct = res_a["pct_of_classified_ipynb"]
            print(
                f"    Done in {res_a['elapsed_seconds']}s | Raw={res_a['counts']['total_raw_reads']:,} | "
                f"Classified={res_a['classified_denominator_ipynb']:,} | "
                f"OnTarget={pa_pct['ontarget_pct']:.2f}% | Nonspec={pa_pct['nonspecific_pct']:.2f}% | "
                f"Dimer={pa_pct['dimer_pct']:.2f}%",
                flush=True,
            )

        print(f"  [Pipeline B] Running standalone Python/edlib on {pid}...", flush=True)
        res_b, mat_b = run_panel_pipeline_b_edlib(
            p, r1_path, r2_path, all_primers, primers_lower, prefix_index,
            amplicons, max_reads=args.max_reads,
        )
        if p["moesm7_not_opt"] is not None:
            res_b["moesm7_matrix_validation"] = compare_with_moesm7(
                mat_b, all_primers, moesm7_path, p["moesm7_not_opt"], p["moesm7_conc_opt"]
            )
            if args.max_reads is None:
                write_pair_table_tsv(
                    mat_b, all_primers, mat_dir / f"{pid.lower()}_pipeline_b_dimer_pairs.tsv.gz",
                    include_directional=True, compress_gz=True,
                )
        elif pid == "E19" and args.max_reads is None:
            write_pair_table_tsv(
                mat_b, all_primers, out_data_dir / "saddle_384_dimer_pairs.tsv",
                include_directional=False, compress_gz=False,
            )

        panel_entry["pipeline_b_edlib"] = res_b
        m1 = res_b["method1_primer_pairing"]
        m2_75 = res_b["method2_edlib_75"]
        m2_80 = res_b["method2_edlib_80"]
        m2_85 = res_b["method2_edlib_85"]
        print(
            f"    Done in {res_b['elapsed_seconds']}s | Classified={res_b['counts']['classified_total']:,} | "
            f"M1={m1['ontarget_pct_classified']:.2f}% | edlib75={m2_75['ontarget_pct_classified']:.2f}% | "
            f"edlib80={m2_80['ontarget_pct_classified']:.2f}% | edlib85={m2_85['ontarget_pct_classified']:.2f}% | "
            f"Dimer={m1['dimer_pct_classified']:.2f}%",
            flush=True,
        )

        canonical["panels"][pid] = panel_entry

    out_json = out_data_dir / "canonical_results.json"
    cc.write_utf8_lf(out_json, json.dumps(canonical, indent=2) + "\n")
    print(f"\nWrote canonical registry -> {out_json}", flush=True)


if __name__ == "__main__":
    main()
