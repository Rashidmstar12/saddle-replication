---
title: "[Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing Validation of Multiplex PCR Read Partitions"
author:
  - name: Rashid Kadayil
    orcid: 0009-0009-6398-4557
    email: rashidkadayil@gmail.com
    affiliations: 1,*
affiliations:
  - code: 1
    name: Independent Researcher, Kerala, India
replication:
  - cite: "Xie, N. G., Wang, M. X., Song, P., Mao, S., Wang, Y., Yang, Y., Luo, J., Ren, S., & Zhang, D. Y. (2022). Designing highly multiplex PCR primer sets with Simulated Annealing Design using Dimer Likelihood Estimation (SADDLE). Nature Communications, 13, 1881."
    bib: "xie2022designing"
    url: "https://doi.org/10.1038/s41467-022-29500-4"
    doi: "10.1038/s41467-022-29500-4"
code:
  - url: "https://github.com/Rashidmstar12/saddle-replication"
  - doi: "10.5281/zenodo.24100087"
data:
  - url: "https://doi.org/10.6084/m9.figshare.16944154.v3"
  - doi: "10.6084/m9.figshare.16944154.v3"
date: 2026-11-10
journal:
  name: "ReScience C"
  issn: 2430-3658
  volume: ""
  issue: ""
  article: ""
---

# Abstract

Multiplex polymerase chain reaction (multiplex PCR) is a foundational enrichment technology for targeted next-generation sequencing (NGS), clinical pathogen panels, and forensic identification. However, high-plex reactions suffer from non-specific primer-dimer formation, which can monopolize sequencing throughput and deplete on-target amplification. In 2022, Xie et al. introduced SADDLE (*Simulated Annealing Design using Dimer Likelihood Estimation*, *Nature Communications*, 13, 1881; PMID: 35410464), describing an empirical heuristic scoring framework and simulated annealing algorithm to design multiplex primer pools up to 384-plex. While the experimental utility of SADDLE was demonstrated in the original paper, the validation sequencing analysis relied on proprietary MATLAB scripts coupled to external command-line aligners, while the simulated annealing design code remained proprietary under patent restrictions. In this work, following the ReScience C guidelines for computational replication, we conduct an independent replication of the central empirical NGS validation finding of Xie et al.: the read classification breakdown and primer-interaction distribution of the SADDLE-designed 384-plex validation panel across 3,942,385 raw paired-end sequencing reads (Sample E19; Figure 3f). Using an independent, open-source Python reimplementation, our pipeline classifies 2,860,658 reads into 43.86% on-target amplicons, 55.27% nonspecific products, and 0.87% primer-dimers, reproducing the published proportions (~43%, ~56%, ~1%) within a 1.0% absolute tolerance. We further corroborate the severe heavy-tailed concentration of primer-dimer formation across the 295,296 candidate primer pairings, where only 1,172 pairs (0.397%) exhibit non-zero dimerization and a single interaction (`rp213` $\times$ `rp371`) accounts for 67.57% of all dimer reads (16,805 reads). In a supplementary exploratory benchmark, we evaluate the predictive validity of the SADDLE Badness heuristic against nearest-neighbor thermodynamics and supervised machine learning models, documenting an Average Precision of 0.0597 (a 4.8-fold improvement over nearest-neighbor free energy alone) while analyzing the severe probability miscalibration of supervised baselines (Brier scores ~0.20 vs ~0.0040 baseline). All replication code, derived datasets, and an automated verification script that executes in under two seconds are openly accessible under permissive licenses.

# Introduction

Multiplex PCR simultaneously amplifies dozens to hundreds of genomic loci within a single reaction volume, serving as an essential preparation step for targeted clinical diagnostics, comprehensive cancer hotspot panels, and metagenomic surveillance [@untergasser2012primer3; @chen2026prism]. However, scaling multiplex PCR beyond small panels is fundamentally constrained by pairwise primer-primer cross-hybridization. For a panel targeting $N$ amplicons, $2N$ distinct primer oligonucleotides are pooled, yielding:

$$\binom{2N+1}{2} = \frac{(2N+1)(2N)}{2} \approx 2N^2$$

possible pairwise interaction combinations (including homodimer and self-interactions). Minor sequence complementarity—specifically when aligned at the extensible 3' terminus—enables DNA polymerase extension, producing short artifactual amplicons termed primer-dimers [@xie2022designing]. Due to their substantially smaller size relative to target genomic amplicons, primer-dimers amplify with high kinetic efficiency, rapidly exhausting enzymatic reagents and saturating sequencing libraries.

To mitigate this exponential risk, Xie et al. [-@xie2022designing] introduced SADDLE (**S**imulated **A**nnealing **D**esign using **D**imer **L**ikelihood **E**stimation), proposing an empirical biophysical scoring heuristic called *Badness* to quantify 3'-directed extensible hybridization and using simulated annealing optimization to partition multiplex primer pools into compatible sub-pools. In their 96-plex experimental optimization (Figures 3a–3c), the authors demonstrated that primer-dimer reads dropped from 90.7% in a naively designed pool (PS1) to 4.9% in an optimized pool (PS3). Scaling the technology to higher multiplexing, the authors designed and experimentally tested a **384-plex primer set designed by SADDLE (768 primers; Sample E19)**. In **Figure 3f** and the accompanying text, the authors reported that the sequencing library of this 384-plex panel achieved high on-target specificity (~43%) and low primer-dimer formation (~1%), with the remaining reads (~56%) representing nonspecific background products aligning across the human genome via non-cognate primer pairs.

Despite the impact of SADDLE, two key reproducibility hurdles exist:
1. **Proprietary Software Barriers:** The simulated annealing software used to design primer sets was withheld under intellectual property and patent restrictions [@xie2022designing]. Furthermore, the empirical NGS validation scripts deposited on Figshare were implemented in proprietary MATLAB and required external aligners [@langmead2012fast], preventing turnkey, open-source execution.
2. **Independent Verification of Empirical Partitions:** While the authors published rounded percentage summaries in Figure 3f of their manuscript, an independent, open-source replication of the read-accounting pipeline operating directly on the raw FASTQ data has not been formally documented.

In this ReScience C article [@rougier2017sustainable], we report an independent computational replication focused strictly on the empirical validation analysis of Xie et al. [-@xie2022designing]:
- **Primary Replication Objective:** Does an independent, open-source Python implementation replicate the published read classification proportions (~43% on-target, ~56% nonspecific, ~1% primer-dimer) across all 3,942,385 raw paired-end reads of the SADDLE 384-plex validation panel (Sample E19, Figure 3f), and does it reproduce the heavy-tailed concentration of primer-dimer reads?
- **Implementation Mapping:** Can the read-classification pipeline be reconstructed without proprietary MATLAB dependencies or external binary alignment tools while maintaining operational fidelity?
- **Supplementary Exploratory Analysis:** In a secondary analysis clearly distinguished from the primary replication claim, how well does the SADDLE Badness heuristic rank observed dimer events compared to standard nearest-neighbor thermodynamics ($-\Delta G$) and supervised machine learning models?

# Original Work and Artifact Audit

The original study by Xie et al. [-@xie2022designing] was published in *Nature Communications* (Volume 13, Article 1881; PMID: 35410464). Primer sequences, genomic coordinates, and assay details were provided in Supplementary Data (MOESM5). The raw paired-end sequencing reads for the 384-plex validation panel (Sample E19, Illumina MiniSeq 2 $\times$ 150 bp) and post-hoc analysis scripts (`Data_process.m`, `Dimer_Identify.m`, `analysis_code.ipynb`) were deposited on Figshare (DOI: 10.6084/m9.figshare.16944154.v3) under a Creative Commons Attribution 4.0 International license (CC BY 4.0).

Our audit of the deposited artifacts identified:
- **Design Code:** The simulated annealing optimization algorithm was not included in the public deposit, with the manuscript stating: *"The code for SADDLE primer design is available upon request for non-commercial academic use"* subject to patent protection. Consequently, the design algorithm itself is out of scope for an open, autonomous replication.
- **Validation Pipeline:** The public analysis workflow comprises MATLAB scripts that trim Illumina adapters, filter reads by Phred quality and length, align reads to reference amplicons using Bowtie2 [@langmead2012fast], match 5' primer sequences via an 8-base hash lookup, and classify reads with insert lengths strictly shorter than $L_1 + L_2 - 2$ bp as primer-dimers.
- **Raw Data Completeness:** The Figshare repository contains the complete 3,942,385 raw paired-end reads for Sample E19 across two gzipped FASTQ files (`E19_384plex_NGS_R1_fastq.gz` and `E19_384plex_NGS_R2_fastq.gz`, totaling 914 MB), providing uncompromised primary empirical inputs.

# Replication Methodology and Implementation

We implemented a clean-room, open-source Python 3.10+ pipeline designed to operate without MATLAB or platform-dependent aligners.

## Structural Code Mapping

Table 1 details the structural correspondence between the original SADDLE validation code and our Python reconstruction.

| Category | Original SADDLE Logic (MATLAB / Bowtie2) | Python Clean-Room Reconstruction | Equivalence Classification |
| :--- | :--- | :--- | :--- |
| **Adapter Trimming** | `Data_process.m:23-57`: Scans `AGATCGGAAGAGCACACGTCTGAACTCC` on R1 and `AGATCGGAAGAGCGTCGTGTAGGGAAAG` on R2; trims read to prefix before adapter. | `reconstruct_saddle_384.py:145-152`: Scans core 17-bp Illumina TruSeq prefixes `AGATCGGAAGAGCACAC` (R1) and `AGATCGGAAGAGCGTCG` (R2); trims reads. | Operationally Equivalent |
| **Quality Filtering** | `Data_process.m:109-111`: Filters read pairs unless $>70\%$ of bases have Phred quality $>30$ (ASCII $>63$), no ambiguous `'N'` bases, and trimmed length $>10$ bp. | `reconstruct_saddle_384.py:155-163`: Exact condition: length $>10$, `'N'` not in sequence, and $\ge 70\%$ bases with quality ASCII $>63$ for both R1 and R2. | Exact Numerical Match |
| **Primer Matching** | `Dimer_Identify.m:28-41, 71-95`: Computes 8-base hash using `PrimeNumberHash(lower(p{i}(1:8)))`; evaluates Smith-Waterman prefix identity $>0.6$. | `reconstruct_saddle_384.py:166-196`: Builds 8-bp prefix index across all 768 primers; evaluates prefix match identity $>0.6$ on prefix of length $\min(\text{len}(p), \text{len}(t))$. | Operationally Equivalent |
| **Insert Length & Dimer Threshold** | `Dimer_Identify.m:96-108`: Dimer condition: `length(currline1) - 1 < (length(f_seq) + length(r_seq) - 2)`. | `reconstruct_saddle_384.py:203-220`: Dimer threshold: $L_1 + L_2 - 2$. Reads with insert length strictly below threshold are classified as dimers. | Exact Formula Match |
| **On-Target Amplicon Assignment** | `analysis_code.ipynb:17-38`: Executes Bowtie2 alignment against reference amplicons; unmapped reads passed to dimer script. | `reconstruct_saddle_384.py:199-210`: Identifies cognate pairs by primer index ($P_2 = P_1 + 384$); classifies cognate reads with insert length $\ge L_1 + L_2 - 2$ as on-target. | Methodological Approximation (Primer pairing vs. sequence alignment) |
| **Pair Aggregation & Homodimers** | `Dimer_Identify.m:98, 125-126`: Increments $768 \times 768$ matrix; preserves homodimers on diagonal `dimer(i, i)`. | `reconstruct_saddle_384.py:238-249`: Computes $768 \times 768$ ordered matrix, then folds into 295,296 unordered pairs with diagonal preserved. | Exact Match |

: Structural mapping and methodological comparison between the published SADDLE validation scripts and the clean-room Python implementation.

## Read Accounting and Partition Logic

Each paired-end sequencing read $(R_1, R_2)$ in Sample E19 is processed through mutually exclusive and exhaustive accounting bins:
1. **Adapter Missing:** Read pairs lacking adapter signatures on both reads, R1 only, or R2 only.
2. **Quality Excluded:** Read pairs with trimmed inserts $\le 10$ bp, containing ambiguous `'N'` bases, or failing the 70% Q30 threshold.
3. **Primer Unmatched:** Read pairs where either or both 5' ends cannot be matched to the 768 pool primers with $>60\%$ identity.
4. **Classified Amplicons:** Reads with both primers successfully resolved, partitioned into **On-Target** (cognate pair with insert $\ge L_1 + L_2 - 2$ bp), **Primer-Dimer** (insert $< L_1 + L_2 - 2$ bp), or **Nonspecific** (non-cognate pairing with insert $\ge L_1 + L_2 - 2$ bp).

# Primary Replication Results

## Replicated NGS Read Accounting Breakdown

Executing the read-accounting pipeline across all 3,942,385 raw paired-end reads of Sample E19 classifies 2,860,658 reads (72.56% of total raw sequencing depth), with 1,081,727 reads (27.44%) excluded. Exclusions partition into:
- Missing adapter sequence: 751,644 reads (19.07% of raw depth);
- Quality and length filters: 137,086 reads (3.48% of raw depth);
- Primer prefix unmatched: 192,997 reads (4.90% of raw depth).

Table 2 compares the classified read proportions directly against the published findings in Figure 3f and main text of Xie et al. [-@xie2022designing].

| Read Category | Replicated Read Count | Replicated Percentage (%) | Published SADDLE Figure 3f [@xie2022designing] | Absolute Difference | Replication Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **On-Target Amplicons** | 1,254,647 | 43.86% | ~43% | +0.86% | **PASS** (within 1.0%) |
| **Nonspecific Background** | 1,581,139 | 55.27% | ~56% | -0.73% | **PASS** (within 1.0%) |
| **Primer-Dimers** | 24,872 | 0.87% | ~1% | -0.13% | **PASS** (within 0.2%) |
| **Total Classified Reads** | 2,860,658 | 100.00% | 100% | 0.00% | **PASS** |

: Comparison between replicated read classification proportions and published SADDLE 384-plex results (Sample E19, MiniSeq 2 $\times$ 150 bp, Figure 3f).

All three categories match the published values within our pre-specified 1.0% absolute tolerance, confirming the quantitative reproducibility of the empirical read partitions.

## Extreme Heavy-Tailed Interaction Concentration

We next analyzed the distribution of the 24,872 observed primer-dimer reads across all 295,296 possible candidate primer pairs:
- **Zero-Dimer Pairs:** 294,124 candidate pairs (99.6031%) produced exactly zero observed dimer reads.
- **Active Dimer Pairs ($\ge 1$ read):** Only 1,172 pairs (0.3969% empirical prevalence) generated one or more dimer reads.
- **Severe Dimer Pairs ($\ge 10$ reads):** Only 98 pairs (0.0332% prevalence) produced 10 or more dimer reads.
- **Dominant Interaction:** A single primer pair—reverse primer 213 (`rp213`) $\times$ reverse primer 371 (`rp371`)—generated 16,805 dimer reads, accounting for **67.57%** of all primer-dimer reads in the entire library.

This reproduces the qualitative finding highlighted by Xie et al.: primer-dimer formation is not uniformly distributed across sequence space, but is characterized by an extreme, heavy-tailed concentration where a tiny minority of catastrophic pairs account for the vast majority of artifactual reads.

# Supplementary Exploratory Analysis: Heuristic and Baseline Benchmarking

*Note: The following benchmark was not published by Xie et al. [-@xie2022designing]. It is conducted here as a secondary exploration to evaluate how the biophysical SADDLE Badness heuristic and machine learning models correlate with empirical dimer labels derived from Sample E19 across all 295,296 candidate primer pairs.*

## Formulation of SADDLE Badness and Baselines

In Xie et al. [-@xie2022designing], Badness is calculated over maximal continuous complementary matches of length $\ge 4$ bp between primers $P_1$ and $P_2$:

$$\text{Badness}(P_1, P_2) = \sum_{k \in \mathcal{M}} 2^{\text{len}_k} \cdot 2^{\text{numGC}_k} \cdot \frac{1}{(d_1^{(k)} + 1)(d_2^{(k)} + 1)}$$

where $\text{len}_k$ is match length, $\text{numGC}_k$ is GC count within the match, and $d_1^{(k)}, d_2^{(k)}$ are nucleotide distances from the 3' termini. In our benchmark, we evaluate an approximate 4–7 bp heuristic variant.

We compare Badness against:
1. **SantaLucia Thermodynamics ($-\Delta G$):** Hybridization free energy calculated at 60 °C using unified nearest-neighbor parameters [@santalucia1998unified].
2. **Logistic Regression (Physics Features):** L2-regularized linear model trained on thermodynamic and alignment features.
3. **XGBoost (Physics Features):** Gradient boosted decision tree ensemble trained on the same feature set.

To ensure deterministic evaluation without tie-breaking instability, candidates are ranked by primary score descending, with ties broken by ascending pair index: $\text{rank\_key}(i) = (-\text{score}_i, \text{pair\_index}_i)$.

## Benchmark Metrics and Calibration Analysis

Table 3 presents the evaluation across all 295,296 candidate pairs against empirical binary labels ($\ge 1$ dimer read).

| Model / Heuristic | Scoring Type | ROC-AUC | PR-AUC (AP) | Precision@10 | Precision@50 | Precision@100 | Brier Score |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **SADDLE Badness** | Biophysical Heuristic | 0.7166 | **0.0597** | **40.0%** | **34.0%** | **24.0%** | N/A |
| **SantaLucia ($-\Delta G$)** | Thermodynamic ($\Delta G$) | 0.7105 | 0.0125 | 20.0% | 10.0% | 6.0% | N/A |
| **Logistic Regression** | Supervised Probabilistic | 0.7317 | 0.0205 | 10.0% | 4.0% | 8.0% | 0.1986 |
| **XGBoost** | Supervised Probabilistic | **0.7564** | 0.0371 | 20.0% | 20.0% | 16.0% | 0.2143 |
| *Constant Baseline* | *Empirical Base Rate* | *0.5000* | *0.0040* | *0.4%* | *0.4%* | *0.4%* | *0.0040* |

: Supplementary benchmark of sequence-based heuristics and supervised models across 295,296 candidate primer pairs from the 384-plex validation panel. Binary labels represent observed dimer formation ($\ge 1$ NGS read). Precision@$K$ evaluated using deterministic secondary key sorting.

### Findings from Supplementary Benchmarking

1. **Heuristic vs. Thermodynamics:** SADDLE Badness achieves an Average Precision of 0.0597, representing a **4.8-fold improvement over nearest-neighbor free energy alone** ($-\Delta G$, AP = 0.0125) and a **4.0-fold improvement in Precision@100** (24.0% vs. 6.0%). This supports the hypothesis of Xie et al. that terminal 3' alignment and heavy GC penalties capture enzymatic extension risk more effectively than global hybridization stability.
2. **Top-Rank Discrimination vs. Supervised Models:** SADDLE Badness outperforms both Logistic Regression (AP = 0.0205, P@100 = 8.0%) and XGBoost (AP = 0.0371, P@100 = 16.0%) in top-rank precision, demonstrating that simple biophysical heuristics remain competitive at the extreme head of the distribution.
3. **Severe Probability Miscalibration of Supervised Models:** Independent recalculation of Brier scores reveals that both Logistic Regression (0.1986) and XGBoost (0.2143) are severely miscalibrated relative to the uninformative baseline ($p(1-p) = 0.003953 \approx 0.0040$). Because the supervised models were trained to maximize ranking separation under class imbalance, their mean predicted probabilities (~0.38) greatly exceed the true empirical prevalence (0.00397), yielding large quadratic calibration penalties.

# Methodological Disclosures

1. **Primer Pairing vs. Bowtie2 Alignment:** Our Python reconstruction classifies on-target amplicons by verifying cognate primer identities ($P_1 = \text{fp}_i, P_2 = \text{rp}_i$) and insert length, omitting Bowtie2 reference genome alignment. While this eliminates external dependencies and reproduces overall read proportions within 1%, it represents an operational approximation of the original pipeline.
2. **Scope Limitation:** This replication focuses strictly on the empirical validation read partitions and heuristic scoring. The proprietary simulated annealing optimization algorithm remains withheld under patent restrictions and is not evaluated.
3. **Distributional Terminology:** We characterize the distribution of primer-dimer read counts as heavy-tailed and highly concentrated; no formal power-law distribution fit (e.g., Clauset-Shalizi-Newman methodology) was conducted.
4. **Data and Code Licensing:** Raw sequencing data from Figshare (DOI: 10.6084/m9.figshare.16944154.v3) are licensed under CC BY 4.0. The replication codebase is distributed under the MIT License.

# Conclusion

We have performed an independent computational replication of the empirical NGS validation results from SADDLE [@xie2022designing]. Using an open-source Python workflow, we confirmed that classified reads from the SADDLE-designed 384-plex validation panel (Figure 3f) partition into 43.86% on-target amplicons, 55.27% nonspecific products, and 0.87% primer-dimers, matching the published values (~43%, ~56%, ~1%) within 1.0% absolute tolerance. We corroborated the severe heavy-tailed concentration of dimer reads, with a single primer pair generating 67.57% of all dimer reads. Finally, our supplementary exploratory benchmark illustrated the predictive utility of the SADDLE Badness heuristic over standard thermodynamics while documenting the calibration characteristics of supervised models.

All data, scripts, and a standalone verification suite are available in the public replication repository.

# Code and Data Availability

- **GitHub Repository:** [https://github.com/Rashidmstar12/saddle-replication](https://github.com/Rashidmstar12/saddle-replication)
- **Zenodo Deposit:** [https://doi.org/10.5281/zenodo.24100087](https://doi.org/10.5281/zenodo.24100087)
- **Primary Data Source:** Figshare (DOI: [10.6084/m9.figshare.16944154.v3](https://doi.org/10.6084/m9.figshare.16944154.v3), CC BY 4.0)

Reviewers can verify all primary and supplementary results in under two seconds:

```bash
git clone https://github.com/Rashidmstar12/saddle-replication.git
cd saddle-replication
python code/verify_replication.py
```

To run end-to-end classification directly from the raw FASTQ sequencing files:
```bash
python code/reproduce_from_fastq.py --full
```

# References
