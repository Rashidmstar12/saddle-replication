#!/usr/bin/env python3
"""
code/build_checksums.py

Generate deterministic UTF-8/LF SHA-256 manifest (CHECKSUMS.tsv) for all
deliverable files in the package (excluding CHECKSUMS.tsv itself, __pycache__,
and temporary LaTeX build artifacts).
"""

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import canonical_classifier as cc

EXCLUDE_DIRS = {"__pycache__", ".git"}
EXCLUDE_EXTS = {".aux", ".bbl", ".bcf", ".blg", ".fdb_latexmk", ".fls", ".out", ".run.xml", ".toc", ".pyc"}
EXCLUDE_NAMES = {"CHECKSUMS.tsv"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    base_dir = Path(__file__).resolve().parent.parent
    rows = []
    for p in sorted(base_dir.rglob("*")):
        if not p.is_file():
            continue
        if any(part in EXCLUDE_DIRS for part in p.parts):
            continue
        rel = p.relative_to(base_dir).as_posix()
        if p.name in EXCLUDE_NAMES or p.suffix in EXCLUDE_EXTS:
            continue
        if p.suffix == ".log" and rel.startswith("article/"):
            continue
        rows.append((rel, sha256_file(p), p.stat().st_size))

    out = base_dir / "CHECKSUMS.tsv"
    lines = ["path\tsha256\tsize_bytes"] + [f"{r}\t{h}\t{s}" for r, h, s in rows]
    cc.write_utf8_lf(out, "\n".join(lines) + "\n")
    print(f"Wrote {len(rows)} entries -> {out}")


if __name__ == "__main__":
    main()
