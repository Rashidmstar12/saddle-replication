# [Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing Validation of Multiplex PCR Read Partitions

**A computational replication of:**  
> Xie, N. G., Wang, M. X., Song, P., Mao, S., Wang, Y., Yang, Y., Luo, J., Ren, S., & Zhang, D. Y. (2022).  
> Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE).  
> *Nature Communications*, 13, 1881.  
> DOI: [10.1038/s41467-022-29500-4](https://doi.org/10.1038/s41467-022-29500-4) | PMID: 35393437

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

## 1. Quick Start for Reviewers (~2 Seconds Verification)

Reviewers can verify the complete computational replication in approximately 1.5 seconds using Python 3.10+ and standard numerical packages (`numpy`, `scikit-learn`):

```bash
# Clone the public repository
git clone https://github.com/Rashidmstar12/saddle-replication.git
cd saddle-replication

# Install dependencies (if not already installed)
pip install numpy scikit-learn

# Run the standalone verifier from the repository root
python code/verify_replication.py
```

### Expected Output
```text
==============================================================================
RESCIENCE C REPLICATION VERIFICATION WORKFLOW
Article: [Re] Computational Replication of SADDLE: Empirical Next-Generation
         Sequencing Read Partitions and Primer Interaction Distribution
         in a 384-Plex PCR Panel
Replication Author: Rashid Kadayil (ORCID: 0009-0009-6398-4557)
Original Study: Xie et al., Nature Communications 13, 1881 (2022)
DOI: 10.1038/s41467-022-29500-4
==============================================================================

[1/4] Verifying Primary Replicated Result: NGS Read Accounting (Sample E19)...
  Total raw paired-end reads: 3,942,385 (expected 3,942,385)
  Total classified reads:     2,860,658 (72.56% of raw depth)
  Total excluded reads:       1,081,727 (27.44% of raw depth)
    - Adapter missing:        751,644 (18.93%)
    - QC filtered (Q30/len):  137,086 (3.48%)
    - Primer unmatched:       192,997 (4.90%)

  Comparison with Published SADDLE Findings (Xie et al. 2022, Figure 4):
    - On-target fraction:   Replicated = 43.86% | Published = ~43% [MATCH within 1.0% tolerance]
    - Nonspecific fraction: Replicated = 55.27% | Published = ~56% [MATCH within 1.0% tolerance]
    - Primer-dimer fraction:Replicated = 0.87%  | Published = ~1%  [MATCH within 0.2% tolerance]
  PASS: Primary Replicated Result (NGS read partition) successfully verified.

[2/4] Verifying Empirical Interaction Concentration across 295,296 Candidate Pairs...
  Total candidate primer pairs:  295,296 (unordered pairs: 768*769/2 = 295,296)
  Sum of observed dimer reads:   24,872 (expected 24,872)
  Nonzero dimer pairs (>= 1):    1,172 (0.3969% prevalence)
  High-confidence pairs (>= 10): 98 (0.0332% prevalence)
  Zero-dimer pairs:              294,124 (99.6031%)
  Dominant pair: row 280576 with 16,805 reads (67.57% of all dimers)
  PASS: Empirical heavy-tailed interaction concentration verified.

[3/4] Verifying Supplementary Exploratory Benchmark (Badness vs Baselines)...
  Empirical Prevalence (mean y_true): 0.003969
  Constant Prevalence Baseline Brier: 0.003953 (~0.0040)
  Logistic Regression Mean Pred Prob: 0.3829 | Brier: 0.198578 [Miscalibrated]
  XGBoost Mean Predicted Probability: 0.3883 | Brier: 0.214317 [Miscalibrated]

  Evaluating deterministic lexsort ranking: primary = -score, secondary = pair_index...
  Model                              ROC-AUC  AP       P@10    P@50    P@100   Brier   Status
  ----------------------------------------------------------------------------------
  SADDLE Badness (Heuristic)           0.7166   0.0597   40.0%   34.0%   24.0%   N/A  [MATCH]
  SantaLucia Thermodynamics (-dG)      0.7105   0.0125   20.0%   10.0%    6.0%   N/A  [MATCH]
  Logistic Regression (Physics)        0.7317   0.0205   10.0%    4.0%    8.0% 0.1986 [MATCH]
  XGBoost (Physics Features)           0.7564   0.0371   20.0%   20.0%   16.0% 0.2143 [MATCH]

  Supplementary Benchmark Observations:
  * SADDLE Badness achieves Average Precision = 0.0597 (4.8x higher than SantaLucia -dG = 0.0125).
  * SADDLE Badness achieves Precision@10 = 40.0% and Precision@100 = 24.0% vs XGBoost (20.0%, 16.0%).
  * Supervised models exhibit severe probability miscalibration (Brier ~0.20 vs ~0.0040 baseline).
  PASS: Supplementary exploratory benchmark verified.

[4/4] Verifying submission member checksums from CHECKSUMS.tsv...
  PASS: All member checksums verified.

==============================================================================
FINAL RESULT: ALL RESCIENCE C REPLICATION CHECKS PASSED (1.42 s)
==============================================================================
```

---

## 2. Targeted Claims & Scope

### Primary Replication Result: 384-plex NGS Read Classification
- **Target:** Figure 4 and main text of Xie et al. (2022).
- **Dataset:** Sample E19, MiniSeq paired-end 2 $\times$ 150 bp, 3,942,385 raw reads (Figshare: `10.6084/m9.figshare.16944154.v3`, CC BY 4.0).
- **Replication Outcome:**
  - On-target amplicons: Published ~43% $\rightarrow$ **Replicated 43.86%** ($\Delta = +0.86\%$, PASS).
  - Nonspecific background: Published ~56% $\rightarrow$ **Replicated 55.27%** ($\Delta = -0.73\%$, PASS).
  - Primer-dimers: Published ~1% $\rightarrow$ **Replicated 0.87%** ($\Delta = -0.13\%$, PASS).
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
│   ├── run_baseline_benchmark.py  # Benchmark evaluation across 295,296 pairs
│   ├── verify_read_accounting.py  # Mutually exclusive read partition audit
│   └── verify_replication.py      # Standalone replication verifier (~1.5s runtime)
└── data/
    ├── baseline_predictions_384.tsv      # 295,296 rows of candidate pairs & predictions
    ├── phase0e_recalculated_metrics.json # Full benchmark metric dictionary
    ├── saddle_384_dimer_pairs.tsv        # Empirical read counts for all 1,172 dimer pairs
    ├── saddle_384_read_accounting.json   # Exhaustive read partition counts
    └── saddle_384_summary_metrics.json   # Summary metrics and comparison tolerances
```

---

## 4. End-to-End Pipeline Execution from Raw FASTQ

The repository includes pre-computed derived datasets so reviewers can verify all results in ~1.5 seconds without downloading the 914 MB raw FASTQ files.

If reviewers wish to execute the complete pipeline directly on raw sequencing data:
1. Download `E19_384plex_NGS_R1_fastq.gz` and `E19_384plex_NGS_R2_fastq.gz` from Figshare ([DOI: 10.6084/m9.figshare.16944154.v3](https://doi.org/10.6084/m9.figshare.16944154.v3)).
2. Download `41467_2022_29500_MOESM5_ESM.xlsx` from the original Nature Communications supplementary information.
3. Place these files in `data/raw/`.
4. Run the live pipeline test:
   ```bash
   python code/verify_replication.py --check-fastq
   ```
5. Or recompute the full 3.94 million read accounting from scratch (~3 minutes runtime):
   ```bash
   python code/verify_read_accounting.py
   ```

---

## 5. Methodological Disclosures

1. **Clean-Room Python Implementation:** This pipeline is an independent Python reconstruction based on the methods described by Xie et al. [-@xie2022designing]. It does not execute the authors' original MATLAB or Bowtie2 scripts.
2. **Proprietary Algorithm Status:** SADDLE's simulated annealing primer selection software remains withheld under patent protection. Only the empirical NGS validation and heuristic scoring components are publicly replicable.
3. **Distributional Characterization:** Primer-dimer formation is described as heavy-tailed and highly concentrated; no formal power-law distribution fit was conducted.
4. **Deterministic Tie-Breaking:** All rankings enforce secondary key sorting by `pair_index` to guarantee exact numerical reproducibility across environments.

---

## 6. License

- Replication manuscript and documentation: [CC-BY 4.0 International](https://creativecommons.org/licenses/by/4.0/).
- Replication code, verifier, and analysis scripts: [MIT License](LICENSE).
- Primary sequencing data from Figshare (DOI: 10.6084/m9.figshare.16944154.v3) are governed by the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
