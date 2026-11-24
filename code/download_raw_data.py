#!/usr/bin/env python3
"""
code/download_raw_data.py

Automated, checksum-verified downloader for all 13 external primary inputs:
  - 8 Illumina MiSeq paired-end FASTQ files from Figshare (CC BY 4.0)
  - 2 Nature Communications Supplementary Data spreadsheets (MOESM5, MOESM7; CC BY 4.0)
  - 3 96-plex reference FASTAs (Ref_1.fasta-Ref_3.fasta) from pinned commit
    1e474abf88f08131544f8cb73dfbc65a5e842cf3 of NinaGXie/SADDLE (retrieved
    programmatically rather than redistributed in this package).

Verifies exact byte size and SHA-256 checksum against data/expected_raw_sha256.json
for every file and exits nonzero on any mismatch.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_url(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    curl_bin = shutil.which("curl") or shutil.which("curl.exe")
    if curl_bin:
        subprocess.run([curl_bin, "-fSL", "--retry", "3", "-o", str(tmp), url], check=True)
    else:
        req = urllib.request.Request(url, headers={"User-Agent": "saddle-replication-downloader/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as out:
            shutil.copyfileobj(resp, out, length=1024 * 1024)
    tmp.replace(dest)


def target_path_for(out_dir: Path, name: str, category: str) -> Path:
    if category == "reference_fasta":
        return out_dir / "reference_fastas" / name
    return out_dir / name


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and SHA-256 verify primary SADDLE inputs.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Destination directory for downloaded raw files")
    parser.add_argument("--refs-only", action="store_true", help="Download/verify only pinned Ref_1-3.fasta files")
    parser.add_argument("--verify-only", action="store_true", help="Verify existing files without downloading")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent.parent
    manifest_path = base_dir / "data" / "expected_raw_sha256.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    files = manifest["files"]
    selected = [
        (name, meta)
        for name, meta in files.items()
        if (not args.refs_only or meta.get("category") == "reference_fasta")
    ]

    ok = True
    for name, meta in selected:
        dest = target_path_for(args.out_dir, name, meta.get("category", "fastq"))
        exp_size = int(meta["size_bytes"])
        exp_sha = meta["sha256"]

        if dest.exists() and dest.stat().st_size == exp_size and sha256_file(dest) == exp_sha:
            print(f"  OK [cached]: {dest} ({exp_size:,} bytes, SHA-256 verified)")
            continue

        if args.verify_only:
            if not dest.exists():
                print(f"  FAIL [missing]: {dest}", file=sys.stderr)
            else:
                print(
                    f"  FAIL [mismatch]: {dest} size={dest.stat().st_size} (exp {exp_size}), "
                    f"sha256={sha256_file(dest)} (exp {exp_sha})",
                    file=sys.stderr,
                )
            ok = False
            continue

        print(f"  Downloading {name} -> {dest} ...", flush=True)
        download_url(meta["url"], dest)
        got_size = dest.stat().st_size
        got_sha = sha256_file(dest)
        if got_size != exp_size or got_sha != exp_sha:
            print(
                f"  FAIL [checksum]: {dest} size={got_size} (exp {exp_size}), sha256={got_sha} (exp {exp_sha})",
                file=sys.stderr,
            )
            ok = False
        else:
            print(f"  OK [downloaded]: {dest} ({got_size:,} bytes, SHA-256 verified)")

    if not ok:
        sys.exit(1)
    print(f"All {len(selected)} selected external inputs verified successfully.")


if __name__ == "__main__":
    main()
