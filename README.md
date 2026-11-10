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
Verifies the complete 295,296-row candidate prediction table, read partitions, label matrices, Brier score calibration, and SHA-256 member checksums without downloading raw FASTQ files. **Requires only Python standard library and standard data science packages (no C-compilers or aligners required)**:

```bash
git clone https://github.com/Rashidmstar12/saddle-replication.git
cd saddle-replication
pip install numpy scikit-learn

python code/verify_replication.py
```

### Modality B: Direct Raw FASTQ Execution
Verifies the classification pipeline directly against raw sequencing reads (Sample E19) using `edlib` and `openpyxl`:

```bash
# Prerequisites for raw FASTQ execution:
pip install openpyxl edlib

# Live sample test on first 50,000 raw reads (~1 second):
python code/reproduce_from_fastq.py --sample 50000

# Complete reconstruction across all 3,942,385 reads (~105 seconds runtime):
python code/reproduce_from_fastq.py --full
```
*(Requires raw FASTQ files `E19_384plex_NGS_R1_fastq.gz` and `E19_384plex_NGS_R2_fastq.gz` and `41467_2022_29500_MOESM5_ESM.xlsx` in `data/raw/`. When raw files are absent, scripts exit with code 2 `SKIPPED_MISSING_INPUT` rather than falsely claiming PASS).*

---

## 2. Replicated Findings & Scope

### Primary Replicated Result: SADDLE 384-plex NGS Read Classification
- **Target:** Figure 3f and main text of Xie et al. (2022).
- **Dataset:** SADDLE-designed 384-plex panel (Sample E19, Illumina MiniSeq 2 $\times$ 150 bp, 3,942,385 raw reads; Figshare: `10.6084/m9.figshare.16944154.v3`, CC BY 4.0).
- **Denominator Reconciliation:**
  - **Derived Classified Read Denominator ($N = 2,860,664$, reconciling with published Figure 3f proportions):**
    *(Note: The exact integer count of classified reads was not reported in Figure 3f; we report both the derived post-filter classified denominator and total raw library output).*
    - Published Figure 3f: **~43% on-target**, **~56% nonspecific**, **~1% primer-dimers** (sum = 100%).
    - Method 1 (Primer Pairing): **43.86% on-target**, **55.27% nonspecific**, **0.87% primer-dimers** ($\Delta \le 0.86\%$, PASS).
    - Method 2 (Gapped Alignment $\ge 75\%$): **43.61% on-target**, **55.52% nonspecific**, **0.87% primer-dimers** ($\Delta \le 0.61\%$, PASS).
    - Method 2 (Gapped Alignment $\ge 80\%$): **42.76% on-target**, **56.37% nonspecific**, **0.87% primer-dimers** ($\Delta \le 0.37\%$, PASS).
  - **Total Raw Read Denominator ($N = 3,942,385$):**
    - On-target: 31.64% to 31.82% ($1,247,496$ to $1,254,647$ reads)
    - Nonspecific: 40.11% to 40.29% ($1,581,139$ to $1,588,290$ reads)
    - Primer-dimers: 0.63% ($24,878$ reads)
    - Excluded reads: 27.44% ($1,081,721$ reads: 19.07% adapter missing, 3.48% low QC/short, 4.90% primer unmatched).
- **Empirical Alignment Concordance & Two-Method Confusion Matrix:**
  - Audited via `code/audit_alignment_and_indels.py` across an initial deterministic sample of the first 5,000 qualifying non-dimer reads from Sample E19 tested against all 384 reference amplicons in both forward and reverse-complement orientations (768 target alignments per read). The sample partitions into:
    - 2,203 cognate reads: 2,193 match own target ($\ge 75\%$), 10 fail due to low sequence quality, and **0 match any different target (0.00% cross-target false positive rate)**.
    - 2,797 non-cognate reads: 126 match target sequences on either strand (4.50% of non-cognates, 2.52% of total sample), while 2,671 fail all 384 targets on both strands (95.50% of non-cognates).
    - Arithmetic Check: $2,193 + 10 + 126 + 2,671 = 5,000$ (exact 100% group partition).
  - **Two-Method Confusion Matrix (Cognate-Anchored Method 2 vs Dual-Strand All-Target Search):**
    - Both On-Target: **2,193 reads** (43.86%)
    - Both Off-Target: **2,681 reads** (53.62%), comprising 2,671 non-cognate reads that fail all targets plus the 10 cognate reads that fail their assigned target and also fail all 384 targets on both strands (10/10 fail both methods; both methods classify them as off-target/nonspecific).
    - Cognate OFF / All-Target ON: **126 reads** (2.52% of sample; non-cognate reads aligning to target amplicons).
    - Cognate ON / All-Target OFF: **0 reads** (0.00%; a read matching its assigned target cannot fail all-target search).
    - **Verified Per-Read Binary Concordance:** $(2,193 + 2,681) / 5,000 = \mathbf{97.48\%}$ (4,874 / 5,000 reads agree).
    - **Verified Per-Read Binary Discordance:** $126 / 5,000 = \mathbf{2.52\%}$.
    - Crucially, this concordance is measured within the 5,000-read sample between cognate-anchored and all-target gapped alignment; it is not a full-dataset confusion matrix and does not claim mathematical equivalence to Bowtie2 reference amplicon mapping.
- **Quantitative Indel Analysis:**
  - CIGAR inspection across 615 audit reads recovered by gapped alignment confirms that **100.0% (615/615) contain indels** (7.80% 1-bp deletion, 0.33% 1-bp insertion, 91.87% multi-indel/stutters; **0.0% substitution-only**). Frameshift indels (e.g. `56=1D9=1X37=`) cause unaligned positional identity to drop to 64.1%, whereas gapped alignment preserves 98.1% identity. All 615 per-read CIGAR alignment records are permanently preserved without truncation in `data/alignment_and_indel_audit.json`. These records confirm that alignment gap operations account for the recovery of reads failing positional comparison, documenting the operational need for gapped alignment rather than proving the biological vs technical cause of each individual indel.
- **Interaction Concentration:** 24,872 to 24,878 observed primer-dimer reads across 295,296 candidate primer pairings. Dimer formation is strongly heavy-tailed: only 1,172 pairs form dimers, and a single dominant interaction (`rp213` $\times$ `rp371`, 16,805 reads) accounts for **67.57%** of all dimer reads in the verified pair-level interaction table ($16,805 / 24,872$).

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
│   ├── audit_alignment_and_indels.py # Empirical dual-strand 384-target concordance & CIGAR audit
│   ├── audit_code_equivalence.py  # Structural comparison between MATLAB and Python logic
│   ├── reconstruct_saddle_384.py  # Full read-classification pipeline from FASTQ
│   ├── reproduce_from_fastq.py    # Turnkey FASTQ streaming pipeline (--sample or --full)
│   ├── run_baseline_benchmark.py  # Benchmark evaluation across 295,296 pairs
│   ├── verify_read_accounting.py  # Mutually exclusive read partition audit
│   └── verify_replication.py      # Standalone replication verifier (~1.5s runtime)
└── data/
    ├── alignment_and_indel_audit.json      # Empirical all-target concordance & 615 CIGAR output
    ├── baseline_predictions_384.tsv        # 295,296 rows of candidate pairs & predictions
    ├── full_fastq_run_audit.log            # Execution audit log across all 3.94M raw FASTQ reads
    ├── phase0e_recalculated_metrics.json   # Full benchmark metric dictionary
    ├── reference_384_amplicons.fasta       # Sequences of all 384 intended target amplicons
    ├── reference_alignment_comparison.json  # Read counts comparing Method 1 vs Method 2
    ├── saddle_384_dimer_pairs.tsv          # Empirical read counts for all 1,172 dimer pairs
    ├── saddle_384_read_accounting.json     # Exhaustive read partition counts
    └── saddle_384_summary_metrics.json     # Summary metrics and comparison tolerances
```

---

## 4. Methodological Disclosures

1. **Replication Scope and Framework:** This work is an independent, open-source **qualified partial replication** and approximate reanalysis of the empirical NGS validation from Xie et al. (2022). While Xie et al. openly provided their post-sequencing analysis workflow (`analysis_code.ipynb`, `Data_process.m`, `Dimer_Identify.m`) and Bowtie2 index files under CC BY 4.0, executing it requires commercial MATLAB and external Bowtie2 installations. In our lightweight, standalone Python implementation, Method 2 executes cognate-anchored semi-global gapped alignment (`edlib`) against the 384 intended target amplicons from MOESM5. Both methods independently confirm the published read proportions (~43% on-target, ~56% nonspecific, ~1% primer-dimer) within 0.86% absolute difference (Method 1 difference = +0.86%, Method 2 difference = +0.61%, both within a 1.0 percentage point tolerance), and dual-strand all-target testing within a 5,000-read audit sample confirms 97.48% operational binary concordance across a verified 2x2 confusion matrix without external dependencies.
2. **Dual Evaluation: Primer Pairing vs. Genuine Gapped Alignment:** This replication provides both a turnkey primer-pairing heuristic (Method 1) and genuine gapped semi-global sequence alignment via `edlib` (Method 2), confirming that both methods replicate published read proportions (~43% on-target, ~56% nonspecific, ~1% primer-dimer) within 0.86% absolute difference (Method 1 difference = +0.86%, Method 2 difference = +0.61%, within a 1.0 percentage point tolerance).
3. **Proprietary Algorithm Status:** SADDLE's simulated annealing primer selection software remains withheld under patent protection. Only the empirical NGS validation and heuristic scoring components are publicly replicable.
4. **Distributional Characterization:** Primer-dimer formation is described as heavy-tailed and highly concentrated; no formal power-law distribution fit was conducted.
5. **Deterministic Tie-Breaking:** All rankings enforce secondary key sorting by `pair_index` to guarantee exact numerical reproducibility across environments.

---

## 5. License

- Replication manuscript and documentation: [CC-BY 4.0 International](https://creativecommons.org/licenses/by/4.0/).
- Replication code, verifier, and analysis scripts: [MIT License](LICENSE).
- Primary sequencing data from Figshare (DOI: 10.6084/m9.figshare.16944154.v3) are governed by the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
