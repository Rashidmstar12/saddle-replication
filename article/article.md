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
bibliography: bibliography.bib
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
4. **Classified Amplicons:** Reads with both primers successfully resolved. These are partitioned into **Primer-Dimers** (insert $< L_1 + L_2 - 2$ bp), and non-dimer products. To evaluate on-target specificity rigorously, non-dimer products are assessed under two distinct methodological criteria:
   - **Method 1 (Primer Pairing Heuristic):** A read is classified as *On-Target* if both primers form a cognate target pair ($P_2 = P_1 + 384$) and insert length $\ge L_1 + L_2 - 2$ bp; all non-cognate pairings are classified as *Nonspecific*. This turnkey method avoids external aligner dependencies.
   - **Method 2 (Genuine Gapped Sequence Alignment via edlib):** In addition to cognate primer matching, the insert sequence is aligned directly against the intended full-length reference amplicon sequence extracted from Supplementary Data MOESM5 [@xie2022designing] using the open-source `edlib` C/Python library (Myers' exact bit-vector semi-global alignment with indels and substitutions, `mode='SHW'`). Reads achieving $\ge 75\%$ or $\ge 80\%$ sequence identity against the intended amplicon are classified as *On-Target*; reads failing this identity threshold are classified as *Nonspecific*.

# Primary Replication Results

## Replicated NGS Read Accounting Breakdown

Executing the read-accounting pipeline across all 3,942,385 raw paired-end reads of Sample E19 classifies 2,860,658 to 2,860,664 reads (72.56% of total raw sequencing depth), with 1,081,721 to 1,081,727 reads (27.44%) excluded. Exclusions partition into:
- Missing adapter sequence: 751,644 reads (19.07% of raw depth);
- Quality and length filters: 137,086 reads (3.48% of raw depth);
- Primer prefix unmatched: 192,991 to 192,997 reads (4.90% of raw depth).

*(Note on the 6-Read Lineage Discrepancy: The minor difference of exactly 6 reads between accounting partitions arises from the denominator used during primer prefix matching. In extreme-short primer dimers with inserts of 12–14 bp, dividing matches by full primer length ($\approx 22$ bp) yields identity $\le 0.60$ and excludes the reads as primer-unmatched [192,997 excluded, 24,872 dimers], whereas dividing by the trimmed insert length evaluates identity as $1.0$ and retains them as dimers [192,991 excluded, 24,878 dimers]. This single-read lineage difference does not affect any published percentage).*

Table 2 compares the classified read proportions directly against the published findings in Figure 3f and the main text of Xie et al. [-@xie2022designing] across all evaluation methodologies and denominators.

| Read Category | Method 1: Primer Pairing Count (% Classified) | Method 2: edlib Gapped ($\ge 75\%$) Count (% Classified) | Method 2: edlib Gapped ($\ge 80\%$) Count (% Classified) | % of Total Raw Reads ($N = 3,942,385$) | Published SADDLE Figure 3f [@xie2022designing] | Replication Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **On-Target Amplicons** | 1,254,647 (**43.86%**) | 1,247,496 (**43.61%**) | 1,223,216 (**42.76%**) | 31.64%–31.82% | ~43% | **PASS** (reproduces ~43%) |
| **Nonspecific Background** | 1,581,139 (**55.27%**) | 1,588,290 (**55.52%**) | 1,612,570 (**56.37%**) | 40.11%–40.29% | ~56% | **PASS** (reproduces ~56%) |
| **Primer-Dimers** | 24,872–24,878 (**0.87%**) | 24,878 (**0.87%**) | 24,878 (**0.87%**) | 0.63% | ~1% | **PASS** (invariant replicate) |
| **Total Classified Reads** | 2,860,658–2,860,664 (**100%**) | 2,860,664 (**100.00%**) | 2,860,664 (**100.00%**) | 72.56% | 100% | **PASS** |
| *Excluded / Filtered Reads* | — | — | — | 27.44% (1,081,721) | — | Excluded by QC / adapter |

: Methodological comparison of replicated read classification proportions against published SADDLE 384-plex results (Sample E19, MiniSeq 2 $\times$ 150 bp, Figure 3f). Method 1 evaluates cognate primer pairing; Method 2 evaluates genuine gapped semi-global sequence alignment via `edlib` against intended reference amplicons. Percentages are presented on both the classified-read basis ($N = 2,860,664$, matching published Figure 3f) and the total raw sequencer output basis ($N = 3,942,385$).

### Biological and Methodological Insights

1. **Resolution of Preprocessing Denominator:** Figure 3f in Xie et al. [-@xie2022designing] presents a partition where the three published categories (~43% on-target, ~56% nonspecific, ~1% primer-dimers) sum to 100%. In next-generation sequencing amplicon analysis, raw sequencer reads that lack adapter sequence (19.07%), fail Phred quality/length filters (3.48%), or fail 5' primer matching (4.90%) are excluded prior to amplicon classification. When evaluated on the classified read denominator ($N = 2,860,664$), our replicated on-target (43.61%–43.86%), nonspecific (55.27%–55.52%), and primer-dimer (0.87%) fractions reproduce the published values within 0.86% absolute difference. When calculated over all raw instrument reads ($N = 3,942,385$), on-target amplicons represent 31.64% to 31.82%, nonspecific products represent 40.11% to 40.29%, primer-dimers represent 0.63%, and excluded reads account for 27.44%. Both denominators are fully documented in `data/reference_alignment_comparison.json` and the execution audit logs.
2. **Invariance of Primer-Dimer Suppression:** Under all classification definitions and denominators, primer-dimer formation accounts for exactly **0.87% of classified reads** (and 0.63% of raw reads), replicating the published finding of **~1%** in Figure 3f within 0.13% absolute difference. Because primer-dimers are defined strictly by insert length truncation ($< L_1 + L_2 - 2$ bp), dimer quantification is completely robust to downstream sequence alignment criteria.
3. **Quantitative Indel Analysis and Frameshift Mechanism:** In earlier exploratory audits, evaluating reference sequence matching via ungapped character-by-character string comparison artificially depressed on-target classification to 39.49% (1,129,696 reads). To determine the exact molecular cause, our audit script (`code/audit_alignment_and_indels.py`) parsed full CIGAR alignment paths across 615 audit reads that passed gapped alignment ($\ge 75\%$ identity) but failed ungapped positional matching. Every single recovered read (**100.0%, 615/615**) contained insertions or deletions: 7.80% (48/615) contained 1-bp deletions, 0.33% (2/615) contained 1-bp insertions, and 91.87% (565/615) exhibited multiple indels or homopolymer stutters; **0.0% had substitutions only**. For example, in read `@M04612:271:000000000-CWC2J:1:1101:16920:1534` matching amplicon 351 (`CIGAR: 56=1D9=1X37=`), a single 1-bp deletion at position 56 caused unaligned positional identity to drop from 98.06% to 64.08% by shifting the downstream reading frame. Gapped semi-global alignment (`edlib`) resolves this frame shift, restoring on-target classification to **43.61%** ($\ge 75\%$ identity, $\Delta = +0.61\%$) and **42.76%** ($\ge 80\%$ identity, $\Delta = -0.24\%$). Complete CIGAR strings, read identifiers, and indel counts are recorded in `data/alignment_and_indel_audit.json`.
4. **Empirical Concordance Between Cognate-Anchored Alignment and All-384-Target Search:** To test whether restricting sequence alignment to the cognate target amplicon introduces classification bias, we developed an independent empirical audit script (`code/audit_alignment_and_indels.py`) and recorded full per-read results in `data/alignment_and_indel_audit.json`. Aligning a representative audit sample of 5,000 classified non-dimer reads against all 384 reference amplicons simultaneously partitions the sample into:
   - *Cognate Reads (2,203 reads, 44.06% of sample):* Exactly 2,193 reads (99.55% of cognates) aligned to their intended cognate target ($\ge 75\%$ identity) and 10 reads (0.45%) failed due to low sequence quality. Crucially, **0 reads (0.00%)** aligned to a different amplicon, demonstrating zero cross-target false positive assignments across the 384 loci.
   - *Non-Cognate Reads (2,797 reads, 55.94% of sample):* Exactly 2,736 reads (97.82% of non-cognates) failed alignment against all 384 amplicons, while 61 reads (2.18% of non-cognates, 1.22% of total sample) achieved $\ge 75\%$ sequence alignment to a target sequence.
   - *Exact Arithmetic Reconciliation:* The sum of all four mutually exclusive partitions ($2,193 + 10 + 61 + 2,736 = 5,000$) reconciles 100% of the tested reads without unexplained records.
   - *Operational Concordance:* Cognate-anchored classification yields 43.86% on-target (2,193 / 5,000), while exhaustive all-384-target alignment yields 45.08% on-target ((2,193 + 61) / 5,000), confirming **98.78% operational concordance** ($100\% - 1.22\%$). Both demonstrate that multiplex PCR pools designed under SADDLE achieve $\approx 43\%$ on-target specificity while rigorously suppressing primer-dimer amplification to under 1%.
5. **Alignment Scope Clarification:** The original study [@xie2022designing] utilized Bowtie2 for genome-wide mapping. In our lightweight, turnkey Python pipeline, Method 2 executes cognate-anchored semi-global gapped alignment (`edlib`, Myers' bit-vector algorithm) against the intended target amplicons from MOESM5. As shown above, this achieves high operational concordance without requiring external binary aligner installations.

## Extreme Heavy-Tailed Interaction Concentration

We next analyzed the distribution of the observed primer-dimer reads across all 295,296 possible candidate primer pairs:
- **Zero-Dimer Pairs:** 294,124 candidate pairs (99.6031%) produced exactly zero observed dimer reads.
- **Active Dimer Pairs ($\ge 1$ read):** Only 1,172 pairs (0.3969% empirical prevalence) generated one or more dimer reads.
- **Severe Dimer Pairs ($\ge 10$ reads):** Only 98 pairs (0.0332% prevalence) produced 10 or more dimer reads.
- **Dominant Interaction:** A single primer pair—reverse primer 213 (`rp213`) $\times$ reverse primer 371 (`rp371`)—generated 16,805 dimer reads, accounting for **67.57%** of dimer reads under the archived 24,872-dimer baseline partition ($16,805 / 24,872$) and **67.55%** under the 24,878-dimer prefix-matched rerun ($16,805 / 24,878$).

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

1. **Dual Evaluation: Primer Pairing vs. Genuine Gapped Sequence Alignment:** Rather than relying solely on primer-pairing heuristics or unaligned string comparisons, our replication evaluated both classification paradigms across all 3,942,385 reads:
   - *Method 1 (Primer Pairing):* Assigns on-target status based on cognate primer IDs ($P_2 = P_1 + 384$) and insert length threshold $\ge L_1 + L_2 - 2$ bp, yielding 43.86% on-target, 55.27% nonspecific, and 0.87% primer-dimers. This reproduces published Figure 3f (~43% / ~56% / ~1%) within 1.0% absolute tolerance.
   - *Method 2 (Genuine Gapped Sequence Alignment via edlib):* Performs Myers' exact bit-vector semi-global sequence alignment against all 384 intended reference amplicons from Supplementary Data MOESM5 [@xie2022designing]. At $\ge 75\%$ identity, reads partition into 43.61% on-target, 55.52% nonspecific, and 0.87% primer-dimers ($\Delta = +0.61\%$). At $\ge 80\%$ identity, reads partition into 42.76% on-target, 56.37% nonspecific, and 0.87% primer-dimers ($\Delta = -0.24\%$). Both aligner cutoffs match published Figure 3f within scientific tolerance.
   - *Positional vs Gapped Alignment Insight:* Naive ungapped string matching artificially depresses on-target reads to 39.49% due to 1–2 bp sequencing indels, whereas genuine gapped alignment robustly tolerates sequencing variations. Crucially, primer-dimer quantification (0.87% vs ~1% published) is completely invariant across all methods.
   - *Alignment Scope and Concordance:* In the original SADDLE publication [@xie2022designing], reads were aligned genome-wide using Bowtie2. In our lightweight, turnkey Python implementation, Method 2 executes cognate-anchored semi-global gapped alignment (`edlib`, Myers' bit-vector algorithm) against the 384 intended target amplicons from MOESM5. An all-target audit across 5,000 classified reads confirmed 0.00% cross-target false positive alignment among cognate reads and a minor 1.22% target capture among non-cognate reads, establishing 98.8% operational concordance between cognate-anchored and all-target alignment without requiring external binary aligners.
2. **Preprocessing Denominators:** As detailed in Section 3, published Figure 3f proportions (~43% / ~56% / ~1%) sum to 100% and correspond to the classified-read denominator ($N = 2,860,664$) after adapter trimming and QC filtering. On the total raw sequencer output basis ($N = 3,942,385$), on-target amplicons represent 31.64% to 31.82%, nonspecific products represent 40.11% to 40.29%, primer-dimers represent 0.63%, and excluded reads account for 27.44%.
3. **Scope Limitation:** This replication focuses strictly on the empirical validation read partitions and heuristic scoring. The proprietary simulated annealing optimization algorithm remains withheld under patent restrictions and is not evaluated.
4. **Distributional Terminology:** We characterize the distribution of primer-dimer read counts as heavy-tailed and highly concentrated; no formal power-law distribution fit (e.g., Clauset-Shalizi-Newman methodology) was conducted.
5. **Data and Code Licensing:** Raw sequencing data from Figshare (DOI: 10.6084/m9.figshare.16944154.v3) are licensed under CC BY 4.0. The replication codebase is distributed under the MIT License.

# Conclusion

We have performed an independent computational replication of the empirical NGS validation results from SADDLE [@xie2022designing]. Using an open-source Python workflow, we confirmed that classified reads from the SADDLE-designed 384-plex validation panel (Figure 3f) replicate published values (~43% on-target, ~56% nonspecific, ~1% primer-dimers). Under primer-pairing classification (Method 1), reads partition into 43.86% on-target, 55.27% nonspecific, and 0.87% primer-dimers. Under genuine gapped sequence alignment against the 384 intended target amplicons (Method 2 via `edlib`), reads partition into 43.61% on-target ($\ge 75\%$ identity) or 42.76% on-target ($\ge 80\%$ identity), 55.52% to 56.37% nonspecific, and 0.87% primer-dimers, reproducing the published ~43% / ~56% / ~1% figures within 0.6% absolute tolerance. Evaluated across all raw sequencing reads, on-target amplicons represent 31.64% of sequencer output. Under all criteria and denominators, primer-dimer formation is strictly suppressed (0.87% of classified reads, 0.63% of raw reads), confirming the core experimental claim of SADDLE. We also corroborated the severe heavy-tailed concentration of dimer reads, with a single primer pair generating 67.57% of dimer reads under the baseline partition (16,805 / 24,872) and 67.55% under the rerun partition (16,805 / 24,878). Finally, our supplementary exploratory benchmark illustrated the predictive utility of the SADDLE Badness heuristic over standard thermodynamics while documenting the calibration characteristics of supervised models.

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
