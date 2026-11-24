# [Re] Open-Source Replication of SADDLE NGS Read Partitions Across 96-Plex and 384-Plex Panels

**Author:** Rashid Kadayil (ORCID: `0009-0009-6398-4557`)  
**Affiliation:** Department of Biochemistry and Molecular Biology, School of Life Sciences, Pondicherry University, Puducherry 605014, India  
**Status:** DRAFT V6-Submission-Ready (author-conducted Kaggle Cloud Linux clean-room reproduction passed; pending investigator evidence-ledger sign-off, final declaration confirmation, and public repository/Zenodo synchronization)  
**Target Journal:** *ReScience C* (synchronized full-content `paper.md` + `paper.bib` AND official `ReScience/template` LaTeX `article/article.pdf` + `article/metadata.yaml`)  
**Replicated Article:** Xie et al. (2022), *Nature Communications* 13, 1881, DOI: `10.1038/s41467-022-29500-4` (PMID: 35410464)  
**Raw Data Source:** Figshare DOI: `10.6084/m9.figshare.16944154.v3` (author: Nina Xie; CC BY 4.0; live API snapshot in `data/figshare_v3_metadata_snapshot.json`)

---

## Summary of Replicated Results

All numeric values below are rendered automatically from `data/canonical_results.json` via `code/render_manuscript.py`.

| Panel | Raw Pairs | Pipeline A (Python + Bowtie2) | Pipeline B (Standalone Python/`edlib` $\ge 75\%$ / $\ge 80\%$) | Published Xie et al. (2022) |
| :--- | ---: | :--- | :--- | :--- |
| **PS1 (96-plex, Fig. 3a)** | 1,823,058 | **90.66%** dimers (2.65% on-target, 6.68% nonspecific) | 90.87% dimers (2.32% on-target) | 90.7% dimers |
| **PS2 (96-plex, Fig. 3a)** | 1,952,269 | **39.62%** dimers (32.78% on-target, 27.60% nonspecific) | 41.34% dimers (29.59% on-target) | 39.6% dimers |
| **PS3 (96-plex, Fig. 3a)** | 1,738,418 | **4.97%** dimers (66.16% on-target, 28.87% nonspecific) | 5.34% dimers (61.71% on-target) | 4.9% dimers |
| **E19 (384-plex, Fig. 3f)** | 3,942,385 | **42.76%** on-target, **56.39%** nonspecific, **0.85%** dimers | $\ge 75\%$: 43.59% / 55.54% / 0.87%<br>$\ge 80\%$: 42.75% / 56.38% / 0.87%<br>$\ge 85\%$: 39.87% / 59.26% / 0.87% | ~43% / ~56% / ~1% |

### Cell-by-Cell Matrix Validation vs. Supplementary Data 3 (`MOESM7`)
All six reconstructed 96-plex directional + unordered pair tables are archived under `data/matrices/`:
- **PS1 (`conc optimized`):** $r = 0.99999$, $\rho = 0.9361$, Active Exact = 65.6%, Jaccard = 0.8765, Norm. L1 = 0.51%, Top-10 = 10/10.
- **PS2 (`conc optimized`):** $r = 0.99999$, $\rho = 0.9048$, Active Exact = 56.7%, Jaccard = 0.8196, Norm. L1 = 0.69%, Top-10 = 10/10.
- **PS3 (`conc optimized`):** $r = 0.99999$, $\rho = 0.8722$, Active Exact = 61.9%, Jaccard = 0.7614, Norm. L1 = 1.02%, Top-10 = 10/10.
- **`not optimized` range (PS1–PS3):** $r = 0.1143$–$0.7251$.

### Matched-Pair Uniform Random Audits on E19 (`seed=20261124`)
Complete 12-column per-read classifications for all 5,000 sampled pairs are archived in `data/matched_pair_audit_5000.tsv.gz` (verified across all 6 confusion matrices, 24/24 cells exact):
- **Paired-end Bowtie2:** **41.26%** (2,063/5,000, identical under both-mapped and proper-pair concordant criteria; 95% CI: 39.90--42.63%).
- **`edlib` $\ge 75\%$:** Cognate = **44.88%** (vs. Bowtie2 concordance 94.54%); all-target = **47.50%** (vs. Bowtie2 concordance 93.76%).
- **`edlib` $\ge 80\%$:** Cognate = **44.08%** (vs. Bowtie2 concordance 95.34%); all-target = **45.44%** (vs. Bowtie2 concordance 95.82%).
- **Indel recovery ($N=20,000$):** 1,856 recovered (9.28%, 95% CI: 8.89--9.69%), 100.0% with alignment gaps (median 8.0 bp, IQR 4.0–9.0 bp).

---

## Reproducible Environment & Automated Data Acquisition

### 1. Recommended Environment Setup (Exact Tested Pins)
1. **Linux / macOS (Conda, recommended):**
   ```bash
   conda env create -f environment.yml
   conda activate saddle-replication
   ```
2. **Exact Python 3.12.4 `pip` lock + system Bowtie2 ($\ge 2.5$):**
   ```bash
   pip install -r requirements-lock-py312.txt
   # Linux: sudo apt-get install bowtie2   |   macOS: brew install bowtie2
   ```
3. **Cross-version Python 3.10–3.14 fallback:**
   ```bash
   pip install -r requirements.txt
   ```
4. **Native Windows note:** Bioconda does not distribute native Windows Bowtie2 binaries; on Windows, either install a native Bowtie2 binary on `PATH` (tested with Bowtie2 2.5.5) or run inside **WSL2 Ubuntu** / **Linux**. Author-conducted clean-room verification transcripts across all 9,456,130 read pairs (performed in separate Linux environments; distinct from third-party verification by an independent researcher) are archived under `data/verification_logs/`:
   - **Author-conducted Kaggle Cloud Linux clean-room reproduction (`Linux 6.18.55+`, fresh `environment.yml`, `Python 3.12.4` / `Bowtie2 2.5.5`, `EXIT_CODE=0`):** `kaggle_cleanroom_full_verify.log`, `kaggle_cleanroom_cross_env_report.json`, and `kaggle_cleanroom_env_freeze.txt`.
   - **Author-conducted WSL2 Ubuntu Linux reproduction (`Linux 6.6.87.2`, `Python 3.14.4` / `Bowtie2 2.5.4`, `EXIT_CODE=0`):** `wsl2_linux_full_verify.log` and `wsl2_linux_cross_env_report.json`.

### 2. Automated Checksum-Verified Download of Primary Inputs (13 Files, ~2.16 GB)
```bash
# Downloads 8 Figshare CC BY 4.0 FASTQs, 2 Nature Communications CC BY 4.0 Excel files,
# and 3 pinned-commit Ref_1-3.fasta references (not redistributed in this ZIP), verifying
# exact byte sizes and SHA-256 checksums against data/expected_raw_sha256.json:
python code/download_raw_data.py --out-dir /path/to/raw
```

### 3. Verification Commands
```bash
# Fast mode: 10/10 synthetic+adversarial unit tests, 4-panel closure, 6 directional 96-plex matrices,
# 2 E19 pair tables, full 6-confusion-matrix (24/24 cell) + unique-ID verification of matched_pair_audit_5000.tsv.gz,
# UTF-8/LF byte-for-byte check of all 4 rendered documents, and SHA-256 manifest:
python code/verify_replication.py

# Full mode: re-executes all 4 panels, 5,000-pair audit, and figures into an isolated temp dir
# (never overwrites data/). Verifies:
#   - 100% exact scientific equivalence on canonical_results.json & random_audit_results.json
#   - 100% exact decompressed SHA-256 match on all 9 pair/audit tables
#   - Exact pixel-identical PNG match on matched Matplotlib/FreeType (or shape-exact + MAE < 0.015
#     & pixel_r > 0.95 visual tolerance across different OS/Matplotlib minor versions)
#   - Non-empty regenerated vector PDF figures
python code/verify_replication.py --full --raw-dir /path/to/raw --threads 4
```

### 4. Compiling the Manuscript (LaTeX or Pandoc Markdown)
```bash
# Official ReScience C LaTeX template -> article/article.pdf
cd article
python yaml-to-latex.py -i metadata.yaml -o metadata.tex
xelatex -interaction=nonstopmode article.tex
biber article
xelatex -interaction=nonstopmode article.tex
xelatex -interaction=nonstopmode article.tex
cd ..

# Synchronized Markdown entrypoint (Pandoc with citation processing, from repository root)
pandoc paper.md --citeproc -o paper_preview.pdf
```

## Licensing & Provenance
See `DATA_AND_LICENSES.md` for full license separation (MIT code, CC BY 4.0 manuscript and Figshare data, and programmatic retrieval of third-party `Ref_1–3.fasta` from pinned upstream commit `1e474abf88f08131544f8cb73dfbc65a5e842cf3`).
