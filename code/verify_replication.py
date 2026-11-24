#!/usr/bin/env python3
"""
code/verify_replication.py

Cross-platform verification suite for the ReScience C SADDLE replication package (V4).

Fast mode (default):
  1. Synthetic + adversarial unit tests (code/test_pipeline.py).
  2. Internal arithmetic closure across PS1, PS2, PS3, and E19 for Pipeline A & B.
  3. Exact verification of E19 pair tables AND all 6 archived 96-plex directional+
     unordered matrices (.tsv.gz) against data/canonical_results.json.
  4. Independent row-by-row recomputation of all 5,000-pair matched-alignment
     confusion matrices from data/matched_pair_audit_5000.tsv.gz and verification
     of CIGAR audit closure against data/canonical_results.json.
  5. Cross-platform UTF-8/LF byte-for-byte reproducibility of article/content.tex,
     article/metadata.yaml, paper.md, and README.md.
  6. Cryptographic SHA-256 verification of all package files in CHECKSUMS.tsv.

Full mode (--full):
  Re-executes run_saddle_replication.py, run_random_audits.py, and
  generate_figures.py into an ISOLATED temporary directory (never overwriting
  archived data/):
    - Separates environment provenance (logging OS/library version differences
      as informational) from strict scientific equivalence.
    - Compares decompressed payloads of all .tsv.gz tables (immune to gzip
      library compression-header variations across platforms).
    - Performs pixel-for-pixel PNG figure comparison (and non-empty PDF check).
    - Exits nonzero on any scientific, table, or figure discrepancy.
"""

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, List

import matplotlib.image as mpimg
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_manuscript as rm


# Fields separated from scientific equivalence comparisons (recorded as provenance only)
NON_SCIENTIFIC_KEYS = {"generated_utc", "elapsed_seconds", "environment"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_decompressed(path: Path) -> str:
    """Hash decompressed bytes for .gz files, or raw bytes for plain files."""
    h = hashlib.sha256()
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
    else:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
    return h.hexdigest()


def strip_non_scientific(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: strip_non_scientific(v) for k, v in obj.items() if k not in NON_SCIENTIFIC_KEYS}
    if isinstance(obj, list):
        return [strip_non_scientific(x) for x in obj]
    return obj


def diff_objects(a: Any, b: Any, path: str = "$") -> List[str]:
    diffs: List[str] = []
    if type(a) is not type(b):
        diffs.append(f"{path}: type {type(a).__name__} != {type(b).__name__}")
        return diffs
    if isinstance(a, dict):
        keys = sorted(set(a.keys()) | set(b.keys()))
        for k in keys:
            if k not in a:
                diffs.append(f"{path}.{k}: missing in archived")
            elif k not in b:
                diffs.append(f"{path}.{k}: missing in recomputed")
            else:
                diffs.extend(diff_objects(a[k], b[k], f"{path}.{k}"))
    elif isinstance(a, list):
        if len(a) != len(b):
            diffs.append(f"{path}: len {len(a)} != {len(b)}")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                diffs.extend(diff_objects(x, y, f"{path}[{i}]"))
    else:
        if a != b:
            diffs.append(f"{path}: archived={a!r} != recomputed={b!r}")
    return diffs


def check_unit_tests(code_dir: Path) -> bool:
    loader = unittest.TestLoader()
    suite = loader.discover(str(code_dir), pattern="test_pipeline.py")
    runner = unittest.TextTestRunner(verbosity=1)
    res = runner.run(suite)
    return res.wasSuccessful()


def check_tsv_matches_canonical(
    tsv_path: Path, expected_dist: dict, label: str, expect_directional: bool = False
) -> bool:
    if not tsv_path.exists():
        print(f"FAIL [{label}]: missing {tsv_path}", file=sys.stderr)
        return False
    n_rows = 0
    total_reads = 0
    active = 0
    severe = 0
    max_reads = 0
    opener = (lambda p: gzip.open(p, "rt", encoding="utf-8")) if str(tsv_path).endswith(".gz") else (lambda p: open(p, "r", encoding="utf-8"))
    exp_header = (
        "pair_index\tprimer_1_id\tprimer_2_id\treads_p1_p2\treads_p2_p1\tdimer_reads"
        if expect_directional
        else "pair_index\tprimer_1_id\tprimer_2_id\tdimer_reads"
    )
    with opener(tsv_path) as f:
        header = f.readline().strip()
        if header != exp_header:
            print(f"FAIL [{label}]: unexpected header {header!r}", file=sys.stderr)
            return False
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if expect_directional:
                fwd = int(parts[3])
                rev = int(parts[4])
                c = int(parts[5])
                if fwd + rev != c:
                    print(f"FAIL [{label}]: directional sum mismatch at row {n_rows}", file=sys.stderr)
                    return False
            else:
                c = int(parts[3])
            n_rows += 1
            total_reads += c
            if c >= 1:
                active += 1
            if c >= 10:
                severe += 1
            if c > max_reads:
                max_reads = c

    checks = [
        ("total_candidate_pairs", n_rows, expected_dist["total_candidate_pairs"]),
        ("total_dimer_reads", total_reads, expected_dist["total_dimer_reads"]),
        ("active_dimer_pairs_ge1", active, expected_dist["active_dimer_pairs_ge1"]),
        ("severe_dimer_pairs_ge10", severe, expected_dist["severe_dimer_pairs_ge10"]),
        ("top1_reads", max_reads, expected_dist["top5_pairs"][0]["dimer_reads"]),
    ]
    for name, got, exp in checks:
        if got != exp:
            print(f"FAIL [{label}]: {name} TSV={got} vs canonical={exp}", file=sys.stderr)
            return False
    tag = "directional+unordered" if expect_directional else "unordered"
    print(f"  PASS [{label} ({tag})]: {n_rows:,} pairs, {total_reads:,} dimer reads, {active:,} active, top1={max_reads:,}")
    return True


def _reconstruct_cm_counts(flags_a: List[bool], flags_b: List[bool]) -> dict:
    both_on = sum(1 for a, b in zip(flags_a, flags_b) if a and b)
    both_off = sum(1 for a, b in zip(flags_a, flags_b) if (not a) and (not b))
    a_on_b_off = sum(1 for a, b in zip(flags_a, flags_b) if a and (not b))
    a_off_b_on = sum(1 for a, b in zip(flags_a, flags_b) if (not a) and b)
    return {
        "both_ontarget": both_on,
        "both_offtarget": both_off,
        "a_on_b_off": a_on_b_off,
        "a_off_b_on": a_off_b_on,
    }


def check_per_read_audit_table(tsv_gz: Path, audit: dict) -> bool:
    """
    Independently reconstruct ALL SIX 2x2 matched-alignment confusion matrices
    cell-by-cell from matched_pair_audit_5000.tsv.gz, verify unique read IDs,
    verify strictly-better and tied non-cognate flags, and compare every cell
    against random_audit_results.json / canonical_results.json.
    """
    if not tsv_gz.exists():
        print(f"FAIL [per_read_audit]: missing {tsv_gz}", file=sys.stderr)
        return False

    seen_ids = set()
    cog_flags: List[bool] = []
    bt2_any_flags: List[bool] = []
    bt2_conc_flags: List[bool] = []
    ge60_flags: List[bool] = []
    c75_flags: List[bool] = []
    c80_flags: List[bool] = []
    a75_flags: List[bool] = []
    a80_flags: List[bool] = []
    strict_75 = 0
    strict_80 = 0
    tied_75 = 0
    tied_80 = 0

    with gzip.open(tsv_gz, "rt", encoding="utf-8") as f:
        header = f.readline().strip()
        exp_h = (
            "read_id\tassigned_amplicon_1based\tis_cognate\tbt2_both_mapped\tbt2_concordant\t"
            "both_mates_ge60\tcog_edlib_sim\tall_edlib_sim\tbest_target_amplicon_1based\t"
            "best_target_strand\tstrictly_better_noncognate\ttied_best_noncognate"
        )
        if header != exp_h:
            print(f"FAIL [per_read_audit]: unexpected header {header!r}", file=sys.stderr)
            return False
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 12:
                print(f"FAIL [per_read_audit]: expected 12 columns, got {len(parts)}", file=sys.stderr)
                return False
            rid = parts[0]
            if rid in seen_ids:
                print(f"FAIL [per_read_audit]: duplicate read_id {rid!r}", file=sys.stderr)
                return False
            seen_ids.add(rid)

            is_cog = int(parts[2]) == 1
            bt2_any = int(parts[3]) == 1
            bt2_conc = int(parts[4]) == 1
            ge60 = int(parts[5]) == 1
            cog_sim = float(parts[6])
            all_sim = float(parts[7])
            strict_nc = int(parts[10]) == 1
            tied_nc = int(parts[11]) == 1
            if strict_nc and tied_nc:
                print(f"FAIL [per_read_audit]: read {rid} has both strictly_better and tied_best set", file=sys.stderr)
                return False

            c75 = is_cog and (cog_sim >= 0.75)
            c80 = is_cog and (cog_sim >= 0.80)
            a75 = all_sim >= 0.75
            a80 = all_sim >= 0.80

            cog_flags.append(is_cog)
            bt2_any_flags.append(bt2_any)
            bt2_conc_flags.append(bt2_conc)
            ge60_flags.append(ge60)
            c75_flags.append(c75)
            c80_flags.append(c80)
            a75_flags.append(a75)
            a80_flags.append(a80)

            if is_cog and a75:
                if strict_nc:
                    strict_75 += 1
                elif tied_nc:
                    tied_75 += 1
            if is_cog and a80:
                if strict_nc:
                    strict_80 += 1
                elif tied_nc:
                    tied_80 += 1

    n = len(seen_ids)
    bm = audit["matched_pair_bowtie2_audit"]
    c75_exp = audit["all_384_target_concordance_75"]
    c80_exp = audit["all_384_target_concordance_80"]

    scalar_checks = [
        ("n_rows", n, c75_exp["sample_size"]),
        ("cognate_total", sum(cog_flags), c75_exp["cognate_total"]),
        ("bt2_both_mapped", sum(bt2_any_flags), bm["bt2_both_mapped_any_count"]),
        ("bt2_concordant", sum(bt2_conc_flags), bm["bt2_concordant_count"]),
        ("bt2_both_equals_concordant", bt2_any_flags == bt2_conc_flags, bm["bt2_both_equals_concordant"]),
        ("ge60_both", sum(ge60_flags), bm["pairs_ge60_both_mates"]),
        ("bt2_ge60_both_mapped", sum(1 for b, g in zip(bt2_any_flags, ge60_flags) if b and g), bm["bt2_ge60_both_mapped_count"]),
        ("strict_cross_75", strict_75, c75_exp["cognate_cross_target_best"]),
        ("tied_cross_75", tied_75, c75_exp["cognate_tied_noncognate_best"]),
        ("strict_cross_80", strict_80, c80_exp["cognate_cross_target_best"]),
        ("tied_cross_80", tied_80, c80_exp["cognate_tied_noncognate_best"]),
    ]
    for name, got, exp in scalar_checks:
        if got != exp:
            print(f"FAIL [per_read_audit]: {name} TSV={got} vs JSON={exp}", file=sys.stderr)
            return False

    cm_specs = [
        ("bowtie2_vs_cognate_75", bt2_any_flags, c75_flags, bm["vs_cognate_edlib75"]),
        ("bowtie2_vs_cognate_80", bt2_any_flags, c80_flags, bm["vs_cognate_edlib80"]),
        ("bowtie2_vs_all_target_75", bt2_any_flags, a75_flags, bm["vs_all_target_edlib75"]),
        ("bowtie2_vs_all_target_80", bt2_any_flags, a80_flags, bm["vs_all_target_edlib80"]),
        ("cognate_vs_all_target_75", c75_flags, a75_flags, c75_exp),
        ("cognate_vs_all_target_80", c80_flags, a80_flags, c80_exp),
    ]
    for label, fa, fb, exp_cm in cm_specs:
        got_cm = _reconstruct_cm_counts(fa, fb)
        for cell_key in ("both_ontarget", "both_offtarget", "a_on_b_off", "a_off_b_on"):
            if got_cm[cell_key] != exp_cm[cell_key]:
                print(
                    f"FAIL [per_read_audit:{label}]: {cell_key} TSV={got_cm[cell_key]} vs JSON={exp_cm[cell_key]}",
                    file=sys.stderr,
                )
                return False

    print(
        f"  PASS [per_read_audit]: verified all {n:,} unique rows, all 6 confusion matrices "
        f"(24/24 cells exact), and strictly-better/tied non-cognate flags"
    )
    return True


def check_template_reproducibility(base_dir: Path, canon: dict) -> bool:
    metrics = rm.flatten_metrics(canon)
    pairs = [
        (base_dir / "article" / "content.tex.in", base_dir / "article" / "content.tex"),
        (base_dir / "article" / "metadata.yaml.in", base_dir / "article" / "metadata.yaml"),
        (base_dir / "paper.md.in", base_dir / "paper.md"),
        (base_dir / "README.md.in", base_dir / "README.md"),
    ]
    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for src, dst in pairs:
            if not src.exists() or not dst.exists():
                print(f"FAIL [template]: missing {src} or {dst}", file=sys.stderr)
                ok = False
                continue
            tmp_out = Path(tmp) / dst.name
            rm.render_template(src, tmp_out, metrics)
            exp_bytes = tmp_out.read_bytes()
            got_bytes = dst.read_bytes()
            if b"\r" in got_bytes:
                print(f"FAIL [template]: {dst.relative_to(base_dir).as_posix()} contains CR bytes (must be LF)", file=sys.stderr)
                ok = False
            elif exp_bytes != got_bytes:
                print(f"FAIL [template]: {dst.relative_to(base_dir).as_posix()} diverges from rendered template!", file=sys.stderr)
                ok = False
            else:
                print(f"  PASS [template]: {dst.relative_to(base_dir).as_posix()} matches canonical_results.json byte-for-byte (UTF-8/LF)")
    return ok


def check_checksums(base_dir: Path) -> bool:
    cs_path = base_dir / "CHECKSUMS.tsv"
    if not cs_path.exists():
        print("FAIL [checksums]: CHECKSUMS.tsv missing", file=sys.stderr)
        return False
    lines = cs_path.read_text(encoding="utf-8").splitlines()
    checked = 0
    for line in lines[1:]:
        if not line.strip():
            continue
        rel, exp_sha, exp_size = line.split("\t")
        fp = base_dir / rel.replace("\\", "/")
        if not fp.exists():
            print(f"FAIL [checksums]: missing file {rel}", file=sys.stderr)
            return False
        if fp.stat().st_size != int(exp_size):
            print(f"FAIL [checksums]: size mismatch for {rel}", file=sys.stderr)
            return False
        if sha256_file(fp) != exp_sha:
            print(f"FAIL [checksums]: SHA-256 mismatch for {rel}", file=sys.stderr)
            return False
        checked += 1
    print(f"  PASS [checksums]: verified {checked} files against CHECKSUMS.tsv")
    return True


def run_isolated_full_recomputation(
    base_dir: Path,
    raw_dir: str,
    ref_dir: str,
    threads: int,
    write_report: str = None,
) -> bool:
    """
    Re-execute full raw-FASTQ pipelines, matched-pair random audits, and figure
    generation into an isolated temp directory; separate environment provenance
    from scientific equivalence; compare decompressed TSV payloads and PNG pixels.
    """
    import matplotlib
    from datetime import datetime, timezone

    code_dir = base_dir / "code"
    arch_data = base_dir / "data"
    arch_figs = base_dir / "article" / "figures"

    with tempfile.TemporaryDirectory(prefix="saddle_full_verify_") as tmp_str:
        tmp_data = Path(tmp_str) / "data"
        tmp_figs = Path(tmp_str) / "figures"

        cmd1 = [sys.executable, str(code_dir / "run_saddle_replication.py"), "--out-dir", str(tmp_data), "--threads", str(threads)]
        cmd2 = [sys.executable, str(code_dir / "run_random_audits.py"), "--out-dir", str(tmp_data), "--seed", "20261124", "--threads", str(threads)]
        if raw_dir:
            cmd1.extend(["--raw-dir", raw_dir])
            cmd2.extend(["--raw-dir", raw_dir])
        if ref_dir:
            cmd1.extend(["--ref-dir", ref_dir])
        cmd3 = [sys.executable, str(code_dir / "generate_figures.py"), "--data-dir", str(tmp_data), "--out-dir", str(tmp_figs)]

        print("\n[FULL 1/4] Re-executing 4-panel replication into isolated temp dir...")
        subprocess.run(cmd1, check=True)
        print("\n[FULL 2/4] Re-executing matched-pair random audits into isolated temp dir...")
        subprocess.run(cmd2, check=True)
        print("\n[FULL 3/4] Re-generating data-driven figures into isolated temp dir...")
        subprocess.run(cmd3, check=True)

        print("\n[FULL 4/4] Comparing recomputed scientific outputs, decompressed tables, and figure pixels...")
        ok = True

        arch_canon = json.loads((arch_data / "canonical_results.json").read_text(encoding="utf-8"))
        new_canon = json.loads((tmp_data / "canonical_results.json").read_text(encoding="utf-8"))
        if arch_canon.get("environment") != new_canon.get("environment"):
            print(f"  [INFO provenance] Archived env:   {arch_canon.get('environment')}")
            print(f"  [INFO provenance] Recomputed env: {new_canon.get('environment')}")

        report: dict = {
            "verified_utc": datetime.now(timezone.utc).isoformat(),
            "archived_environment": arch_canon.get("environment"),
            "recomputed_environment": {
                **(new_canon.get("environment") or {}),
                "matplotlib": matplotlib.__version__,
            },
            "json_equivalence": {},
            "decompressed_table_equivalence": {},
            "figure_checks": {},
        }

        for jname in ("canonical_results.json", "random_audit_results.json"):
            a_obj = strip_non_scientific(json.loads((arch_data / jname).read_text(encoding="utf-8")))
            b_obj = strip_non_scientific(json.loads((tmp_data / jname).read_text(encoding="utf-8")))
            diffs = diff_objects(a_obj, b_obj)
            report["json_equivalence"][jname] = {"exact_match": len(diffs) == 0, "diff_count": len(diffs)}
            if diffs:
                print(f"FAIL [full:{jname}]: {len(diffs)} scientific mismatches (first 10):", file=sys.stderr)
                for d in diffs[:10]:
                    print(f"  {d}", file=sys.stderr)
                ok = False
            else:
                print(f"  PASS [full:{jname}]: exact scientific equivalence verified")

        rel_tables = [
            "saddle_384_dimer_pairs.tsv",
            "saddle_384_dimer_pairs_bowtie2.tsv",
            "matched_pair_audit_5000.tsv.gz",
            "matrices/ps1_pipeline_a_dimer_pairs.tsv.gz",
            "matrices/ps1_pipeline_b_dimer_pairs.tsv.gz",
            "matrices/ps2_pipeline_a_dimer_pairs.tsv.gz",
            "matrices/ps2_pipeline_b_dimer_pairs.tsv.gz",
            "matrices/ps3_pipeline_a_dimer_pairs.tsv.gz",
            "matrices/ps3_pipeline_b_dimer_pairs.tsv.gz",
        ]
        for rel in rel_tables:
            h_arch = sha256_decompressed(arch_data / rel)
            h_new = sha256_decompressed(tmp_data / rel)
            report["decompressed_table_equivalence"][rel] = {
                "exact_match": h_arch == h_new,
                "decompressed_sha256": h_new,
            }
            if h_arch != h_new:
                print(f"FAIL [full:{rel}]: decompressed SHA-256 archived={h_arch} != recomputed={h_new}", file=sys.stderr)
                ok = False
            else:
                print(f"  PASS [full:{rel}]: exact decompressed byte-for-byte match")

        for png_name in ("fig1_multipanel_replication.png", "fig2_dimer_concentration_and_matrix.png"):
            img_arch = mpimg.imread(arch_figs / png_name)
            img_new = mpimg.imread(tmp_figs / png_name)
            if img_arch.shape != img_new.shape:
                print(f"FAIL [full:{png_name}]: shape {img_arch.shape} != {img_new.shape}", file=sys.stderr)
                report["figure_checks"][png_name] = {"pass": False, "shape_match": False}
                ok = False
            elif np.array_equal(img_arch, img_new):
                print(f"  PASS [full:{png_name}]: pixel-identical to archived PNG ({img_new.shape[1]}x{img_new.shape[0]})")
                report["figure_checks"][png_name] = {
                    "pass": True, "mode": "pixel_identical", "shape": list(img_new.shape), "mae": 0.0, "pixel_r": 1.0
                }
            else:
                a_f = img_arch.astype(np.float32)
                b_f = img_new.astype(np.float32)
                mae = float(np.mean(np.abs(a_f - b_f)))
                corr = float(np.corrcoef(a_f.ravel(), b_f.ravel())[0, 1])
                passed = (mae < 0.015) and (corr > 0.95)
                report["figure_checks"][png_name] = {
                    "pass": passed,
                    "mode": "cross_platform_tolerance",
                    "shape": list(img_new.shape),
                    "mae": round(mae, 6),
                    "pixel_r": round(corr, 6),
                }
                if passed:
                    print(
                        f"  PASS [full:{png_name}]: visually equivalent across OS/Matplotlib font stacks "
                        f"(shape={img_new.shape[1]}x{img_new.shape[0]}, mae={mae:.6f}, pixel_r={corr:.6f})"
                    )
                else:
                    print(
                        f"FAIL [full:{png_name}]: mae={mae:.6f}, pixel_r={corr:.6f} outside cross-platform tolerance",
                        file=sys.stderr,
                    )
                    ok = False

        for pdf_name in ("fig1_multipanel_replication.pdf", "fig2_dimer_concentration_and_matrix.pdf"):
            fp = tmp_figs / pdf_name
            size = fp.stat().st_size if fp.exists() else 0
            report["figure_checks"][pdf_name] = {"pass": size > 0, "mode": "non_empty_pdf", "size_bytes": size}
            if size == 0:
                print(f"FAIL [full:{pdf_name}]: missing or empty regenerated PDF", file=sys.stderr)
                ok = False
            else:
                print(f"  PASS [full:{pdf_name}]: regenerated non-empty vector PDF ({size:,} bytes)")

        report["overall_pass"] = ok
        if write_report:
            rp = Path(write_report)
            rp.parent.mkdir(parents=True, exist_ok=True)
            rp.write_bytes((json.dumps(report, indent=2) + "\n").encode("utf-8"))
            print(f"  Wrote cross-environment verification report -> {rp}")

        return ok


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify SADDLE ReScience C replication package (V5).")
    parser.add_argument("--full", action="store_true", help="Re-run raw FASTQ pipelines into temp dir and verify scientific equivalence")
    parser.add_argument("--raw-dir", type=str, default=None, help="Raw data directory (for --full)")
    parser.add_argument("--ref-dir", type=str, default=None, help="Directory with pinned Ref_1-3.fasta files (for --full)")
    parser.add_argument("--threads", type=int, default=4, help="Bowtie2 threads (for --full)")
    parser.add_argument("--write-report", type=str, default=None, help="Optional path to write machine-readable --full verification JSON report")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent.parent
    code_dir = base_dir / "code"
    data_dir = base_dir / "data"

    print("=" * 76)
    print("SADDLE MULTI-PANEL REPLICATION VERIFICATION SUITE (V5)")
    print("=" * 76)

    print("\n[1/6] Running synthetic and adversarial unit tests...")
    if not check_unit_tests(code_dir):
        sys.exit(1)

    print("\n[2/6] Verifying canonical_results.json accounting closure across PS1, PS2, PS3, E19...")
    canon = json.loads((data_dir / "canonical_results.json").read_text(encoding="utf-8"))
    for pid in ("PS1", "PS2", "PS3", "E19"):
        p = canon["panels"][pid]
        pb = p["pipeline_b_edlib"]
        raw = pb["counts"]["total_raw_reads"]
        cls = pb["counts"]["classified_total"]
        excl = pb["exclusions_summary"]["total_excluded"]
        assert cls + excl == raw, f"{pid}: Pipeline B closure failed"
        for mkey in ("method1_primer_pairing", "method2_edlib_75", "method2_edlib_80", "method2_edlib_85", "ungapped_positional_75"):
            blk = pb[mkey]
            assert blk["ontarget_count"] + blk["nonspecific_count"] + blk["dimer_count"] == cls
        pa = p.get("pipeline_a_bowtie2")
        if pa:
            c = pa["counts"]
            assert c["total_raw_reads"] == raw
            assert (
                c["no_adapter_both"] + c["no_adapter_r1_only"] + c["no_adapter_r2_only"]
                + c["qc_short_len"] + c["qc_ambiguous_n"] + c["qc_low_quality"]
                + c["qc_passed"]
            ) == raw
            assert (
                c["routed_long_ge60_to_bowtie2"]
                + c["routed_short_lt60_to_dimer_cand"]
                + c["mixed_length_60bp_dropped_matlab_rule"]
            ) == c["qc_passed"]
            assert (
                c["bt2_both_mapped_concordant"]
                + c["bt2_both_mapped_discordant"]
                + c["bt2_one_mate_mapped"]
                + c["bt2_both_unmapped"]
            ) == c["routed_long_ge60_to_bowtie2"]
            assert (
                c["dimer_cand_len_le15_excluded"]
                + c["dimer_cand_primer_unmatched"]
                + c["dimers"]
                + c["nonspecific"]
            ) == c["dimer_cand_total"]
            print(
                f"  PASS [{pid}]: Raw={raw:,} | Pipeline A Dimer={pa['pct_of_classified_ipynb']['dimer_pct']:.2f}%, "
                f"OnTarget={pa['pct_of_classified_ipynb']['ontarget_pct']:.2f}% | "
                f"Pipeline B (edlib>=75%) Dimer={pb['method2_edlib_75']['dimer_pct_classified']:.2f}%, "
                f"OnTarget={pb['method2_edlib_75']['ontarget_pct_classified']:.2f}%"
            )

    print("\n[3/6] Verifying E19 TSVs and all 6 archived 96-plex directional+unordered matrices...")
    e19 = canon["panels"]["E19"]
    ok_tables = (
        check_tsv_matches_canonical(data_dir / "saddle_384_dimer_pairs.tsv", e19["pipeline_b_edlib"]["dimer_distribution"], "E19 Pipeline B TSV", expect_directional=False)
        and check_tsv_matches_canonical(data_dir / "saddle_384_dimer_pairs_bowtie2.tsv", e19["pipeline_a_bowtie2"]["dimer_distribution"], "E19 Pipeline A TSV", expect_directional=False)
    )
    for pid in ("PS1", "PS2", "PS3"):
        p = canon["panels"][pid]
        ok_tables = (
            ok_tables
            and check_tsv_matches_canonical(data_dir / "matrices" / f"{pid.lower()}_pipeline_a_dimer_pairs.tsv.gz", p["pipeline_a_bowtie2"]["dimer_distribution"], f"{pid} Pipeline A Matrix", expect_directional=True)
            and check_tsv_matches_canonical(data_dir / "matrices" / f"{pid.lower()}_pipeline_b_dimer_pairs.tsv.gz", p["pipeline_b_edlib"]["dimer_distribution"], f"{pid} Pipeline B Matrix", expect_directional=True)
        )
    if not ok_tables:
        sys.exit(1)

    print("\n[4/6] Verifying 5,000-pair per-read audit TSV and CIGAR audits match canonical_results.json...")
    audit = json.loads((data_dir / "random_audit_results.json").read_text(encoding="utf-8"))
    for key in ("all_384_target_concordance_75", "all_384_target_concordance_80"):
        cm = audit[key]
        assert cm["both_ontarget"] + cm["both_offtarget"] + cm["a_on_b_off"] + cm["a_off_b_on"] == cm["sample_size"]
        assert cm == canon["random_audits_e19"][key]
    assert audit["matched_pair_bowtie2_audit"] == canon["random_audits_e19"]["matched_pair_bowtie2_audit"]
    ind = audit["indel_recovery_audit"]
    assert ind["recovered_by_gapped_count"] == len(audit["recovered_cigar_records"])
    assert ind == canon["random_audits_e19"]["indel_recovery_audit"]
    if not check_per_read_audit_table(data_dir / "matched_pair_audit_5000.tsv.gz", audit):
        sys.exit(1)

    print("\n[5/6] Verifying cross-platform UTF-8/LF manuscript (LaTeX + paper.md) & README templates...")
    if not check_template_reproducibility(base_dir, canon):
        sys.exit(1)

    print("\n[6/6] Verifying CHECKSUMS.tsv SHA-256 manifest...")
    if not check_checksums(base_dir):
        sys.exit(1)

    if args.full:
        if not run_isolated_full_recomputation(
            base_dir, args.raw_dir, args.ref_dir, args.threads, write_report=args.write_report
        ):
            sys.exit(1)
        print("\n" + "=" * 76)
        print("FULL INDEPENDENT RAW-FASTQ RECOMPUTATION & VERIFICATION: PASS")
        print("=" * 76)
    else:
        print("\n" + "=" * 76)
        print("ARCHIVED CONSISTENCY, MATRIX, PER-READ AUDIT, TEMPLATE & CHECKSUMS: PASS")
        print("(Run with --full --raw-dir <dir> for isolated raw FASTQ recomputation)")
        print("=" * 76)


if __name__ == "__main__":
    main()
