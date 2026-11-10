# [Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing Validation of Multiplex PCR Read Partitions

**A computational replication of:**  
> Xie, N. G., Wang, M. X., Song, P., Mao, S., Wang, Y., Yang, Y., Luo, J., Ren, S., & Zhang, D. Y. (2022).  
> Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE).  
> *Nature Communications*, 13, 1881 (Figure 3f).  
> DOI: [10.1038/s41467-022-29500-4](https://doi.org/10.1038/s41467-022-29500-4) | PMID: 35410464

**Replication Author:**  
Rashid Kadayil ([ORCID: 0009-0009-6398-4557](https://orcid.org/0009-0009-6398-4557))  
Independent Researcher, Kerala, India  
Email: `rashidkadayil@gmail.com`

**Repository & Identifiers:**
- **Public GitHub Repository:** [https://github.com/Rashidmstar12/saddle-replication](https://github.com/Rashidmstar12/saddle-replication)
- **Zenodo Archive:** [10.5281/zenodo.24100087](https://doi.org/10.5281/zenodo.24100087) (Concept DOI: [10.5281/zenodo.24100086](https://doi.org/10.5281/zenodo.24100086))
- **Primary Data Source:** Figshare [10.6084/m9.figshare.16944154.v3](https://doi.org/10.6084/m9.figshare.16944154.v3) (CC BY 4.0)
- **ReScience C Metadata:** `article/metadata.yaml`

---

## 1. Quick Start for Reviewers

Reviewers can verify this computational replication through two independent modalities:

### Modality A: Fast Standalone Consistency Check (~1.5 Seconds)
Verifies the complete 295,296-row candidate prediction table, read partitions, label matrices, Brier score calibration, and SHA-256 member checksums without downloading raw FASTQ files:

```bash
git clone https://github.com/Rashidmstar12/saddle-replication.git
cd saddle-replication
pip install numpy scikit-learn openpyxl edlib

python code/verify_replication.py
```

### Modality B: Direct Raw FASTQ Execution
Verifies the classification pipeline directly against raw sequencing reads (Sample E19):

```bash
# Live sample test on first 50,000 raw reads (~1 second):
python code/reproduce_from_fastq.py --sample 50000

# Complete reconstruction across all 3,942,385 reads (~2 minutes runtime):
python code/reproduce_from_fastq.py --full
```
*(Requires raw FASTQ files `E19_384plex_NGS_R1_fastq.gz` and `E19_384plex_NGS_R2_fastq.gz` and `41467_2022_29500_MOESM5_ESM.xlsx` in `data/raw/`. When raw files are absent, scripts exit with code 2 `SKIPPED_MISSING_INPUT` rather than falsely claiming PASS).*

---

## 2. Replicated Findings & Scope

### Primary Replicated Result: SADDLE 384-plex NGS Read Classification
- **Target:** Figure 3f and main text of Xie et al. (2022).
- **Dataset:** SADDLE-designed 384-plex panel (Sample E19, Illumina MiniSeq 2 $\times$ 150 bp, 3,942,385 raw reads; Figshare: `10.6084/m9.figshare.16944154.v3`, CC BY 4.0).
- **Replication Outcomes Across Both Methodologies:**
  - **Method 1 (Primer Pairing):**
    - On-target amplicons: Published Figure 3f ~43% $\rightarrow$ **Replicated 43.86%** ($\Delta = +0.86\%$, PASS).
    - Nonspecific background: Published Figure 3f ~56% $\rightarrow$ **Replicated 55.27%** ($\Delta = -0.73\%$, PASS).
    - Primer-dimers: Published Figure 3f ~1% $\rightarrow$ **Replicated 0.87%** ($\Delta = -0.13\%$, PASS).
  - **Method 2 (Genuine Gapped Sequence Alignment via edlib):**
    - At $\ge 75\%$ sequence identity: **43.61% on-target** ($1,247,496$ reads), **55.52% nonspecific** ($1,588,290$ reads), **0.87% dimers** ($\Delta = +0.61\%$, PASS).
    - At $\ge 80\%$ sequence identity: **42.76% on-target** ($1,223,216$ reads), **56.37% nonspecific** ($1,612,570$ reads), **0.87% dimers** ($\Delta = -0.24\%$, PASS).
    - Primer-dimer suppression is **100% invariant** across all methods at **0.87%** (matching published ~1%).
  - **Diagnostic Lineage Notes:**
    - *Ungapped String Comparison:* Ungapped positional comparison artificially depressed on-target classification to 39.49% due to 1–2 bp sequencing indels, whereas genuine gapped alignment (`edlib`) resolves these variations.
    - *6-Read Difference:* Exactly 6 extreme-short primer dimer reads (12–14 bp) were excluded by full-length primer denominator (`matches/len(p) <= 0.60`) in earlier scripts, but retained as dimers when evaluated over trimmed insert prefixes.
- **Interaction Concentration:** 24,872 observed primer-dimer reads across 295,296 candidate primer pairings. Dimer formation is strongly heavy-tailed: only 1,172 pairs form dimers, and a single dominant interaction (`rp213` $\times$ `rp371`) accounts for 16,805 reads (**67.57%** of all primer-dimer reads).

### Supplementary Exploratory Benchmark: Heuristics vs. Baselines
*Note: This benchmark was not published in Xie et al. (2022); it is presented as a supplementary exploratory analysis.*
- SADDLE Badness provides a **4.8-fold improvement in Average Precision** (0.0597 vs. 0.0125) and a **4.0-fold improvement in Precision@100** (24.0% vs. 6.0%) over nearest-neighbor free energy alone ($-\Delta G$).
- Supervised models (Logistic Regression, XGBoost) were optimized for ranking discrimination under severe class imbalance (0.397% prevalence), resulting in high Brier scores (0.1986, 0.2143) due to probability miscalibration (mean prediction ~0.38 vs. base rate 0.0040).

---

## 3. Directory Layout

```text
saddle-replication/
├── CHECKSUMS.tsv                  # Cryptographic SHA-256 manifest of all files
├── LICENSE                        # MIT License
├── README.md                      # This reviewer guide
├── article/
│   ├── article.html               # Formatted standalone HTML with citations
│   ├── article.md                 # Manuscript source in Markdown (ReScience format)
│   ├── article.tex                # Generated LaTeX source
│   ├── bibliography.bib           # Verified BibTeX references with authentic DOIs/PMIDs
│   └── metadata.yaml              # Official ReScience C metadata file
├── code/
│   ├── audit_code_equivalence.py  # Structural comparison between MATLAB and Python logic
│   ├── reconstruct_saddle_384.py  # Full read-classification pipeline from FASTQ
│   ├── reproduce_from_fastq.py    # Turnkey FASTQ streaming pipeline (--sample or --full)
│   ├── run_baseline_benchmark.py  # Benchmark evaluation across 295,296 pairs
│   ├── verify_read_accounting.py  # Mutually exclusive read partition audit
│   └── verify_replication.py      # Standalone replication verifier (~1.5s runtime)
└── data/
    ├── baseline_predictions_384.tsv        # 295,296 rows of candidate pairs & predictions
    ├── phase0e_recalculated_metrics.json   # Full benchmark metric dictionary
    ├── reference_384_amplicons.fasta       # Sequences of all 384 intended target amplicons
    ├── reference_alignment_comparison.json  # Read counts comparing Method 1 vs Method 2
    ├── saddle_384_dimer_pairs.tsv          # Empirical read counts for all 1,172 dimer pairs
    ├── saddle_384_read_accounting.json     # Exhaustive read partition counts
    └── saddle_384_summary_metrics.json     # Summary metrics and comparison tolerances
```

---

## 4. Methodological Disclosures

1. **Dual Evaluation: Primer Pairing vs. Genuine Gapped Alignment:** This replication provides both a turnkey primer-pairing heuristic (Method 1) and genuine gapped semi-global sequence alignment via `edlib` (Method 2), confirming that both methods replicate published read proportions (~43% on-target, ~56% nonspecific, ~1% primer-dimer) within 0.6% absolute tolerance.
2. **Proprietary Algorithm Status:** SADDLE's simulated annealing primer selection software remains withheld under patent protection. Only the empirical NGS validation and heuristic scoring components are publicly replicable.
3. **Distributional Characterization:** Primer-dimer formation is described as heavy-tailed and highly concentrated; no formal power-law distribution fit was conducted.
4. **Deterministic Tie-Breaking:** All rankings enforce secondary key sorting by `pair_index` to guarantee exact numerical reproducibility across environments.

---

## 5. License

- Replication manuscript and documentation: [CC-BY 4.0 International](https://creativecommons.org/licenses/by/4.0/).
- Replication code, verifier, and analysis scripts: [MIT License](LICENSE).
- Primary sequencing data from Figshare (DOI: 10.6084/m9.figshare.16944154.v3) are governed by the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
