---
title: "[Re] Computational Replication of SADDLE: Empirical Next-Generation Sequencing Read Partitions and Primer Interaction Distribution in a 384-Plex PCR Panel"
author:
  - Rashid Kadayil
authors:
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

Multiplex polymerase chain reaction (multiplex PCR) is a foundational enrichment technology for targeted next-generation sequencing (NGS), clinical pathogen panels, and forensic identification. However, high-plex reactions suffer from non-specific primer-dimer formation, which can monopolize sequencing throughput and deplete on-target amplification. In 2022, Xie et al. introduced SADDLE (*Simulated Annealing Design using Dimer Likelihood Estimation*, *Nature Communications*, 13, 1881; PMID: 35410464), describing an empirical heuristic scoring framework and simulated annealing algorithm to design multiplex primer pools up to 384-plex. While the original authors made their post-sequencing analysis scripts (`analysis_code.ipynb`, `Data_process.m`, `Dimer_Identify.m`) publicly available through GitHub and deposited raw sequencing data on Figshare under a CC BY 4.0 license, executing their workflow requires commercial MATLAB and external Bowtie2 installations, while the proprietary simulated annealing primer design algorithm was withheld under patent protection. In this work, following ReScience C guidelines, we present an independent **qualified partial replication** focused on the central empirical NGS validation finding of Xie et al.: the read classification breakdown and primer-interaction distribution of the SADDLE-designed 384-plex validation panel across 3,942,385 raw paired-end sequencing reads (Sample E19; Figure 3f). Using an independent, open-source Python reimplementation incorporating both primer-pairing classification (Method 1) and genuine gapped semi-global sequence alignment against assigned target amplicons (Method 2 via `edlib`), our pipeline classifies a derived post-filter denominator of 2,860,664 reads into 43.61% to 43.86% on-target amplicons, 55.27% to 55.52% nonspecific products, and 0.87% primer-dimers (24,878 reads), approximately reproducing the published Figure 3f proportions (~43%, ~56%, ~1%) within 0.86% for Method 1 and 0.61% for Method 2 (within a 1.0 percentage point scientific tolerance). Evaluated over the full raw instrument library ($N = 3,942,385$), on-target amplicons represent 31.64% of sequencer output. We further corroborate the severe heavy-tailed concentration of primer-dimer formation across the 295,296 candidate primer pairings, where only 1,172 pairs (0.397%) exhibit non-zero dimerization and a single dominant interaction (`rp213` $\times$ `rp371`) accounts for 67.57% of all dimer reads in the verified pair-level interaction matrix (16,805 / 24,872 reads). In a supplementary exploratory benchmark, we evaluate the predictive validity of the SADDLE Badness heuristic against nearest-neighbor thermodynamics and supervised machine learning models, documenting an Average Precision of 0.0597 (a 4.8-fold improvement over nearest-neighbor free energy alone) while analyzing the severe probability miscalibration of supervised baselines (Brier scores ~0.20 vs ~0.0040 baseline). All replication code, derived datasets, and an automated verification script that executes in under two seconds are openly accessible under permissive licenses.

# Introduction

Multiplex PCR simultaneously amplifies dozens to hundreds of genomic loci within a single reaction volume, serving as an essential preparation step for targeted clinical diagnostics, comprehensive cancer hotspot panels, and metagenomic surveillance [@untergasser2012primer3; @chen2026prism]. However, scaling multiplex PCR beyond small panels is fundamentally constrained by pairwise primer-primer cross-hybridization. For a panel targeting $N$ amplicons, $2N$ distinct primer oligonucleotides are pooled, yielding:

$$\binom{2N+1}{2} = \frac{(2N+1)(2N)}{2} \approx 2N^2$$

possible pairwise interaction combinations (including homodimer and self-interactions). Minor sequence complementarity—specifically when aligned at the extensible 3' terminus—enables DNA polymerase extension, producing short artifactual amplicons termed primer-dimers [@xie2022designing]. Due to their substantially smaller size relative to target genomic amplicons, primer-dimers amplify with high kinetic efficiency, rapidly exhausting enzymatic reagents and saturating sequencing libraries.

To mitigate this quadratically increasing interaction burden, Xie et al. [-@xie2022designing] introduced SADDLE (**S**imulated **A**nnealing **D**esign using **D**imer **L**ikelihood **E**stimation), proposing an empirical biophysical scoring heuristic called *Badness* to quantify 3'-directed extensible hybridization and using simulated annealing optimization to partition multiplex primer pools into compatible sub-pools. In their 96-plex experimental optimization (Figures 3a–3c), the authors demonstrated that primer-dimer reads dropped from 90.7% in a naively designed pool (PS1) to 4.9% in an optimized pool (PS3). Scaling the technology to higher multiplexing, the authors designed and experimentally tested a **384-plex primer set designed by SADDLE (768 primers; Sample E19)**. In **Figure 3f** and the accompanying text, the authors reported that the sequencing library of this 384-plex panel achieved high on-target specificity (~43%) and low primer-dimer formation (~1%), with the remaining reads (~56%) representing nonspecific background products formed by non-cognate primer pairings.

Despite the impact of SADDLE, computational reproducibility is constrained by two key factors:
1. **Software Ecosystem and Proprietary Restrictions:** While the original authors made their post-sequencing analysis workflow publicly available on GitHub and Figshare (comprising a Jupyter notebook `analysis_code.ipynb`, MATLAB scripts `Data_process.m` and `Dimer_Identify.m`, and Bowtie2 index files), with raw sequencing data deposited under CC BY 4.0, executing this workflow requires a commercial MATLAB license and an external Bowtie2 installation configured on Linux/macOS. Furthermore, SADDLE's simulated annealing primer pool design software was withheld under intellectual property and patent protection ("available upon request for non-commercial academic use" [@xie2022designing]).
2. **Independent Verification of Empirical Partitions:** While the authors published rounded percentage summaries in Figure 3f of their manuscript (~43% on-target, ~56% nonspecific, ~1% primer-dimers), an independent, open-source replication of the read-accounting pipeline operating directly on the raw FASTQ data in a standalone, turnkey environment has not been formally documented.

In this ReScience C article [@rougier2017sustainable], we report an independent **qualified partial replication and alternative Python reanalysis** focused on the empirical validation analysis of Xie et al. [-@xie2022designing]:
- **Primary Replication Objective:** To assess whether an independent, open-source Python reanalysis of the publicly available 384-plex SADDLE sequencing data (Sample E19) yields read-classification proportions consistent with those reported by Xie et al. [-@xie2022designing, Figure 3f], while explicitly documenting preprocessing exclusions, read-count denominators, and methodological differences from the original MATLAB/Bowtie2 workflow.
- **Secondary Objective:** To characterize the distribution and concentration of observed primer-dimer interactions within the reconstructed read-count data.
- **Supplementary Exploratory Analysis:** To compare approximate SADDLE-style primer-dimer risk scores with thermodynamic and supervised machine-learning baselines, without presenting these comparisons as replicated findings from the original study.

# Original Work and Artifact Audit

The original study by Xie et al. [-@xie2022designing] was published in *Nature Communications* (Volume 13, Article 1881; PMID: 35410464). Primer sequences, genomic coordinates, and assay details were provided in Supplementary Data (MOESM5). The raw paired-end sequencing reads for the 384-plex validation panel (Sample E19, Illumina MiSeq paired-end sequencing using MiSeq Reagent Kit v2) were deposited on Figshare (DOI: 10.6084/m9.figshare.16944154.v3) under a Creative Commons Attribution 4.0 International license (CC BY 4.0), while analysis scripts (`Data_process.m`, `Dimer_Identify.m`, `analysis_code.ipynb`) were made publicly available through GitHub and Figshare.

Our audit of the deposited artifacts identified:
- **Design Code:** The simulated annealing optimization algorithm was not included in the public deposit, with the manuscript stating: *"The code for SADDLE primer design is available upon request for non-commercial academic use"* subject to patent protection. Consequently, the design algorithm itself is out of scope for an open, autonomous replication.
- **Validation Pipeline:** The authors made their post-hoc analysis scripts (`analysis_code.ipynb`, `Data_process.m`, `Dimer_Identify.m`) and Bowtie2 reference indices publicly available on GitHub and Figshare. The pipeline trims Illumina adapters, filters reads by Phred quality and length, aligns reads to intended reference amplicons using Bowtie2 [@langmead2012fast], matches 5' primer sequences via an 8-base hash lookup, and classifies reads with insert lengths strictly shorter than $L_1 + L_2 - 2$ bp as primer-dimers.
- **Raw Data Completeness:** The Figshare repository contains the complete 3,942,385 raw paired-end reads for Sample E19 across two gzipped FASTQ files (`E19_384plex_NGS_R1_fastq.gz` and `E19_384plex_NGS_R2_fastq.gz`, totaling 914 MB), providing uncompromised primary empirical inputs.

# Replication Methodology and Implementation

We implemented a clean-room, open-source Python 3.10+ pipeline designed to operate without MATLAB or platform-dependent aligners.

## Structural Code Mapping

Table 1 details the structural correspondence and methodological differences between the published SADDLE validation code (`Data_process.m`, `Dimer_Identify.m`, Bowtie2) and our Python reconstruction.

| Category | Original SADDLE Logic (MATLAB / Bowtie2) | Python Clean-Room Reconstruction | Methodological Relationship |
| :--- | :--- | :--- | :--- |
| **Adapter Trimming & Routing** | `Data_process.m:23-57, 118, 128`: Scans full 28–29 bp adapters (`AGATCGGAAGAGCACACGTCTGAACTCC` on R1, `AGATCGGAAGAGCGTCGTGTAGGGAAAG` on R2); bifurcates reads by length ($\ge 60$ bp routed to Bowtie2 vs. $< 60$ bp routed to candidate dimer search). | `reconstruct_saddle_384.py:145-152`: Scans core 17-bp Illumina TruSeq prefixes (`AGATCGGAAGAGCACAC` on R1, `AGATCGGAAGAGCGTCG` on R2); processes all qualifying reads in a single unified streaming flow without 60-bp length bifurcation. | **Alternative Python Heuristic** (17-bp prefix matching; single-stream flow without 60-bp length bifurcation) |
| **Quality Filtering** | `Data_process.m:109-111`: Filters read pairs unless $>70\%$ of bases have Phred quality $>30$ (ASCII $>63$), no ambiguous `'N'` bases, and trimmed length $>10$ bp. | `reconstruct_saddle_384.py:155-163`: Exact condition: length $>10$ bp, `'N'` not in sequence, and $> 70\%$ bases with quality ASCII $>63$ for both R1 and R2. | **Threshold Reimplementation** ($>70\%$ Q30, no 'N', length $>10$ bp) |
| **Primer Matching** | `Dimer_Identify.m:28-41, 71-95`: 8-base hash lookup via `PrimeNumberHash`; computes Smith–Waterman local alignment score via `swalign` (`score / score_perfect > 0.6`). | `reconstruct_saddle_384.py:166-196`: 8-bp prefix index across all 768 primers; evaluates direct positional match fraction $>0.60$ on prefix of length $\min(\text{len}(p), \text{len}(t))$. | **Alternative Heuristic** (Positional character matching vs. MATLAB `swalign` Smith–Waterman local alignment) |
| **Insert Length & Dimer Threshold** | `Dimer_Identify.m:70, 96-108`: Requires `length(readF) > 15` and `length(f_seq) < length(currline1)`; dimer condition: insert length $< L_1 + L_2 - 2$ bp. | `reconstruct_saddle_384.py:203-220`: Dimer condition: insert length $< L_1 + L_2 - 2$ bp; does not apply MATLAB's $> 15$ bp filter. | **Formula Approximation** (Shared threshold $< L_1 + L_2 - 2$ bp; differs in read length $> 15$ bp eligibility filter) |
| **On-Target Amplicon Assignment** | `analysis_code.ipynb:17-38`: Executes Bowtie2 alignment against reference amplicons; unmapped reads routed to dimer script. | `reconstruct_saddle_384.py:199-210`: Method 1 assigns cognate pairs by primer ID ($P_2 = P_1 + 384$); Method 2 executes cognate-anchored semi-global gapped alignment (`edlib`) against assigned reference amplicon. | **Methodological Approximation** (Cognate primer pairing and edlib semi-global gapped alignment vs. Bowtie2) |
| **Pair Aggregation & Panel Scaling** | `Dimer_Identify.m:3, 52, 98, 125-126`: Public repository script hardcoded for 192 primers (`index = 1:192; % 96-plex`), populating $192 \times 192$ matrix. | `reconstruct_saddle_384.py:238-249`: Constructs $768 \times 768$ interaction matrix across all 768 primers for the 384-plex panel, then folds into 295,296 unordered pairs. | **Panel Scaling Reimplementation** ($768 \times 768$ matrix for 384-plex; public MATLAB script hard-coded for $192 \times 192$ / 96-plex) |

: Structural mapping and methodological comparison between the published SADDLE validation scripts and the clean-room Python implementation.

## Read Accounting and Partition Logic

Each paired-end sequencing read $(R_1, R_2)$ in Sample E19 is processed through mutually exclusive and exhaustive accounting bins:
1. **Adapter Missing:** Read pairs lacking adapter signatures on both reads, R1 only, or R2 only.
2. **Quality Excluded:** Read pairs with trimmed inserts $\le 10$ bp, containing ambiguous `'N'` bases, or failing the 70% Q30 threshold.
3. **Primer Unmatched:** Read pairs where either or both 5' ends cannot be matched to the 768 pool primers with $>60\%$ identity.
4. **Classified Amplicons:** Reads with both primers successfully resolved. These are partitioned into **Primer-Dimers** (insert $< L_1 + L_2 - 2$ bp), and non-dimer products. To evaluate on-target specificity rigorously, non-dimer products are assessed under two distinct methodological criteria:
   - **Method 1 (Primer Pairing Heuristic):** A read is classified as *On-Target* if both primers form a cognate target pair ($P_2 = P_1 + 384$) and insert length $\ge L_1 + L_2 - 2$ bp; all non-cognate pairings are classified as *Nonspecific*. This turnkey method avoids external aligner dependencies.
   - **Method 2 (Cognate-Anchored Gapped Sequence Alignment via edlib):** In addition to cognate primer matching, the insert sequence is aligned directly against the intended full-length reference amplicon sequence extracted from Supplementary Data MOESM5 [@xie2022designing] using the open-source `edlib` C/Python library [@sosic2017edlib] (Myers' exact bit-vector semi-global alignment with indels and substitutions, `mode='SHW'`). Reads achieving semi-global edit distance similarity $\ge 75\%$ or $\ge 80\%$ (defined as $1 - d_E / L_{\text{read}} \ge 0.75$, where $d_E$ is Myers' bit-vector edit distance and $L_{\text{read}}$ is the trimmed read sequence length, `len(t1)`) against their assigned cognate amplicon are classified as *On-Target*; reads failing this threshold are classified as *Nonspecific*.

# Primary Replication Results

## Replicated NGS Read Accounting Breakdown

Executing the read-accounting pipeline across all 3,942,385 raw paired-end reads of Sample E19 classifies 2,860,658 to 2,860,664 reads (72.56% of total raw sequencing depth), with 1,081,721 to 1,081,727 reads (27.44%) excluded. Exclusions partition into:
- Missing adapter sequence: 751,644 reads (19.07% of raw depth);
- Quality and length filters: 137,086 reads (3.48% of raw depth);
- Primer prefix unmatched: 192,991 to 192,997 reads (4.90% of raw depth).

*(Note on the 6-Read Lineage Discrepancy and Length Eligibility Rules: The minor difference of exactly 6 reads between accounting partitions arises from the denominator used during primer prefix matching in our Python pipeline. In extreme-short primer dimers with inserts of 12–14 bp, dividing matches by full primer length ($\approx 22$ bp) yields identity $\le 0.60$ and excludes the reads as primer-unmatched [192,997 excluded, 24,872 dimers], whereas dividing by the trimmed insert length evaluates identity as $1.0$ and retains them as dimers [192,991 excluded, 24,878 dimers]. Crucially, the original published MATLAB script `Dimer_Identify.m` enforces an explicit length filter requiring read length $> 15$ bp (`if length(readF) > 15`) and insert longer than primer length. Consequently, in the original MATLAB pipeline, reads of 12–14 bp would be filtered out by this length condition. These 6 reads therefore reflect our Python pipeline's alternative length filter ($> 10$ bp) rather than exact reproduction of the original MATLAB eligibility rules. Neither accounting alters the rounded 0.87% dimer fraction).*

Table 2 compares the classified read proportions directly against the published findings in Figure 3f and the main text of Xie et al. [-@xie2022designing] across all evaluation methodologies and denominators.

| Read Category | Method 1: Primer Pairing Count (% Classified) | Method 2: edlib Gapped ($\ge 75\%$) Count (% Classified) | Method 2: edlib Gapped ($\ge 80\%$) Count (% Classified) | % of Total Raw Reads ($N = 3,942,385$) | Published SADDLE Figure 3f [@xie2022designing] | Replication Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **On-Target Amplicons** | 1,254,647 (**43.86%**) | 1,247,496 (**43.61%**) | 1,223,216 (**42.76%**) | 31.64%–31.82% | ~43% | **PASS** (reproduces ~43%) |
| **Nonspecific Background** | 1,581,139 (**55.27%**) | 1,588,290 (**55.52%**) | 1,612,570 (**56.37%**) | 40.11%–40.29% | ~56% | **PASS** (reproduces ~56%) |
| **Primer-Dimers** | 24,872–24,878 (**0.87%**) | 24,878 (**0.87%**) | 24,878 (**0.87%**) | 0.63% | ~1% | **PASS** (invariant replicate) |
| **Total Classified Reads** | 2,860,658–2,860,664 (**100%**) | 2,860,664 (**100.00%**) | 2,860,664 (**100.00%**) | 72.56% | 100% | **PASS** |
| *Excluded / Filtered Reads* | — | — | — | 27.44% (1,081,721) | — | Excluded by QC / adapter |

: Methodological comparison of replicated read classification proportions against published SADDLE 384-plex results (Sample E19, Illumina MiSeq paired-end sequencing with MiSeq Reagent Kit v2, Figure 3f). Method 1 evaluates cognate primer pairing; Method 2 evaluates cognate-anchored gapped semi-global sequence alignment via `edlib` against assigned cognate reference amplicons. Percentages are presented on both the derived post-filter classified-read basis ($N = 2,860,664$, which yields proportions approximately reconciling with published Figure 3f) and the total raw sequencer output basis ($N = 3,942,385$). The exact integer denominator was not published in Figure 3f and remains unverified.

### Biological and Methodological Insights

1. **Derived Preprocessing Denominator and Read Partitioning:** In Figure 3f of Xie et al. [-@xie2022designing], published proportions (~43% on-target, ~56% nonspecific, ~1% primer-dimers) sum to 100%, but the original publication did not disclose the underlying integer read count or filtering denominator. We report two transparent denominators: (a) our pipeline's derived post-filter classified denominator ($N = 2,860,664$ read pairs), after filtering 27.44% unclassifiable reads (19.07% missing adapters, 3.48% low quality/short, 4.90% unmatched primers), which yields on-target (43.61% for Method 2, 43.86% for Method 1), nonspecific (55.27%–55.52%), and primer-dimer (0.87%) fractions in approximate agreement with the published values within a 1.0 percentage point scientific tolerance (differences of +0.61% and +0.86% on-target); and (b) the total raw sequencer output ($N = 3,942,385$ read pairs), across which on-target amplicons represent 31.64% to 31.82%, nonspecific products represent 40.11% to 40.29%, primer-dimers represent 0.63%, and excluded reads represent 27.44%. We do not present the original denominator as established, but rather report this derived post-filter denominator as an approximate operational alignment under our filtering criteria. Both denominators are fully documented in `data/reference_alignment_comparison.json` and the execution audit logs.
2. **Invariance of Primer-Dimer Suppression:** Under all classification definitions and denominators, primer-dimer formation accounts for exactly **0.87% of classified reads** (and 0.63% of raw reads), replicating the published finding of **~1%** in Figure 3f within 0.13% absolute difference. Because primer-dimers are defined strictly by insert length truncation ($< L_1 + L_2 - 2$ bp), dimer quantification is completely robust to downstream sequence alignment criteria.
3. **Quantitative Indel Analysis and Operational Gap Mechanism:** In earlier exploratory audits, evaluating reference sequence matching via ungapped character-by-character string comparison artificially depressed on-target classification to 39.49% (1,129,696 reads). To determine why gapped alignment recovers reads that fail positional matching, our audit script (`code/audit_alignment_and_indels.py`) parsed full CIGAR alignment paths across 615 audit reads that passed gapped semi-global alignment ($\ge 75\%$ similarity) but failed ungapped positional matching. Every single recovered read (**100.0%, 615/615**) contained alignment gaps (insertions or deletions relative to the reference; **0.0% had substitutions only**). Across these 615 records, 53 reads (8.62%) contained $\le 2$ inserted/deleted bases (48 with 1-bp deletion, 2 with 1-bp insertion, 3 with 2-bp indels), while the majority (562 reads, 91.38%) exhibited larger or multiple indels, with a median of 8 inserted/deleted bases across the recovered set. For example, in read `@M04612:271:000000000-CWC2J:1:1101:16920:1534` matching amplicon 351 (`CIGAR: 56=1D9=1X37=`), a 1-bp deletion at position 56 caused unaligned positional identity to drop from 98.06% to 64.08% by shifting the downstream alignment frame. Gapped semi-global alignment (`edlib`) resolves this frame shift, restoring on-target classification to **43.61%** ($\ge 75\%$ similarity, $\Delta = +0.61\%$) and **42.76%** ($\ge 80\%$ similarity, $\Delta = -0.24\%$). Complete CIGAR strings, read identifiers, and indel counts for all 615 recovered reads are recorded without truncation in `data/alignment_and_indel_audit.json`. These alignment records confirm that alignment gap operations account for the recovery of reads failing positional comparison. However, this sample audit reflects a specific subset of reads and does not establish that small 1–2 bp indels or homopolymer stutters are the predominant cause of read divergence across the full 3.94-million-read library, nor does it establish the biological versus technical etiology of individual alignment gaps.
4. **Dual-Strand All-384-Target Alignment Audit and Concordance:** In the full FASTQ reconstruction pipeline (Method 2), each cognate-primer read is aligned semi-globally against its single assigned target amplicon. To test whether this cognate anchoring introduces classification bias compared to an exhaustive search, we developed an independent empirical audit script (`code/audit_alignment_and_indels.py`) and recorded summary confusion statistics, aggregate partition counts, and all 615 recovered indel CIGAR records in `data/alignment_and_indel_audit.json`. Testing an initial deterministic sample of the first 5,000 qualifying non-dimer reads from Sample E19 against all 384 reference amplicons in both forward and reverse-complement orientations (768 target alignments per read) partitions the sample into:
   - *Cognate Reads (2,203 reads, 44.06% of sample):* Exactly 2,193 reads (99.55% of cognates) aligned to their intended cognate target ($\ge 75\%$ identity) and 10 reads (0.45%) failed the 75% alignment threshold against their assigned target. Crucially, **0 reads (0.00%)** had a non-cognate amplicon as their highest-identity alignment $\ge 75\%$ (0/2,203 cross-target best matches; secondary sub-threshold alignments were not logged).
   - *Non-Cognate Reads (2,797 reads, 55.94% of sample):* Exactly 2,671 reads (95.50% of non-cognates) failed alignment against all 384 amplicons on both strands, while 126 reads (4.50% of non-cognates, 2.52% of total sample) achieved $\ge 75\%$ sequence alignment to a target amplicon on either strand.
   - *Exact Arithmetic Reconciliation:* The sum of all four mutually exclusive partitions ($2,193 + 10 + 126 + 2,671 = 5,000$) reconciles 100% of the tested reads without unexplained records.
   - *Two-Method Confusion Matrix and Verified Concordance:* Cross-classifying each read's binary outcome under Method 2 (cognate-anchored gapped alignment) against dual-strand all-384-target alignment constructs a complete 2 $\times$ 2 confusion matrix:
     - *Both On-Target (Agreement):* Exactly 2,193 reads (43.86% of sample).
     - *Both Off-Target (Agreement):* Exactly 2,681 reads (53.62% of sample), comprising 2,671 non-cognate reads failing all targets plus all 10 cognate reads that failed alignment to their assigned target. A dedicated per-read trace confirmed that 100% (10/10) of these cognate reads also failed alignment against all 384 targets on both strands (0 passed any target), meaning both methods agree they are off-target/nonspecific.
     - *Cognate OFF / All-Target ON (Discordant):* Exactly 126 reads (2.52% of sample; non-cognate reads achieving sequence alignment to a target amplicon).
     - *Cognate ON / All-Target OFF (Discordant):* Exactly 0 reads (0.00%; mathematically impossible as the assigned target is within the 384-target set).
     - *Verified Per-Read Binary Concordance:* Total agreement is $(2,193 + 2,681) / 5,000 = \mathbf{97.48\%}$ ($4,874 / 5,000$ reads), with exactly 126 discordant reads ($\mathbf{2.52\%}$). Crucially, this concordance metric measures agreement specifically within this 5,000-read sample between cognate-anchored and all-target semi-global gapped alignment; it is not derived from a full-dataset confusion matrix, should not be extrapolated as proof of identical classification across all 3.94 million reads, and does not establish mathematical equivalence to the original Bowtie2 amplicon mapping workflow.
5. **Replication Scope and Framework Framing:** This study is an independent, open-source **qualified partial replication and alternative Python reanalysis**. In the original publication [@xie2022designing], reads were processed using MATLAB scripts (`Data_process.m`, `Dimer_Identify.m`) and Bowtie2. In our standalone Python workflow, we reconstruct the core classification logic using open-source tools: core adapter prefix trimming, single-stream processing, direct prefix matching, and cognate-anchored semi-global gapped alignment (`edlib`). While our implementation reproduces the reported read classification proportions (~43% / ~56% / ~1%) and the extreme heavy-tailed primer-dimer concentration (~67.6% from one pair) without commercial software licenses or external aligner dependencies, it is an alternative approximate workflow rather than an exact reimplementation of the original MATLAB scripts. SADDLE's simulated annealing primer pool design software was withheld under patent protection and is not evaluated.

## Extreme Heavy-Tailed Interaction Concentration

We next analyzed the distribution of the observed primer-dimer reads across all 295,296 possible candidate primer pairs:
- **Zero-Dimer Pairs:** 294,124 candidate pairs (99.6031%) produced exactly zero observed dimer reads.
- **Active Dimer Pairs ($\ge 1$ read):** Only 1,172 pairs (0.3969% empirical prevalence) generated one or more dimer reads.
- **Severe Dimer Pairs ($\ge 10$ reads):** Only 98 pairs (0.0332% prevalence) produced 10 or more dimer reads.
- **Dominant Interaction:** A single primer pair—reverse primer 213 (`rp213`) $\times$ reverse primer 371 (`rp371`)—generated 16,805 dimer reads, accounting for **67.57%** of all observed dimer reads in the verified pair-level distribution ($16,805 / 24,872$ reads in `data/saddle_384_dimer_pairs.tsv`). In the current rerun partition ($24,878$ total dimer reads), the 6 additional extreme-short reads across the library do not materially alter this extreme concentration (~67.6%).

This reproduces the qualitative finding highlighted by Xie et al.: primer-dimer formation is not uniformly distributed across sequence space, but is characterized by an extreme, heavy-tailed concentration where a tiny minority of catastrophic pairs account for the vast majority of artifactual reads.

# Supplementary Exploratory Analysis: Heuristic and Baseline Benchmarking

*Note: The following benchmark was not published by Xie et al. [-@xie2022designing]. It is conducted here as a secondary exploration to evaluate how the biophysical SADDLE Badness heuristic and machine learning models correlate with empirical dimer labels derived from Sample E19 across all 295,296 candidate primer pairs.*

## Formulation of SADDLE Badness and Baselines

In Xie et al. [-@xie2022designing], Badness is calculated over maximal continuous complementary matches of length $\ge 4$ bp between primers $P_1$ and $P_2$:

$$\text{Badness}(P_1, P_2) = \sum_{k \in \mathcal{M}} 2^{\text{len}_k} \cdot 2^{\text{numGC}_k} \cdot \frac{1}{(d_1^{(k)} + 1)(d_2^{(k)} + 1)}$$

where $\text{len}_k$ is match length, $\text{numGC}_k$ is GC count within the match, and $d_1^{(k)}, d_2^{(k)}$ are nucleotide distances from the 3' termini. In our benchmark, we evaluate an approximate 4–7 bp heuristic variant.

We compare Badness against:
1. **SantaLucia Thermodynamics ($-\Delta G$):** Approximate nearest-neighbor hybridization free-energy score based on standard 37 °C parameters ($\Delta G_{37}^\circ$) [@santalucia1998unified].
2. **Logistic Regression (Augmented Features):** L2-regularized linear model trained on thermodynamic free energy ($\Delta G_{37}^\circ$), terminal 3' hybridization, match lengths, GC content, and the approximate SADDLE-style Badness heuristic score.
3. **XGBoost (Augmented Features):** Gradient boosted decision tree ensemble trained on the exact same augmented feature set (thermodynamics, alignment features, and approximate SADDLE-style Badness).

To ensure deterministic evaluation without tie-breaking instability, candidates are ranked by primary score descending, with ties broken by ascending pair index: $\text{rank\_key}(i) = (-\text{score}_i, \text{pair\_index}_i)$.

## Benchmark Metrics and Calibration Analysis

Table 3 presents the evaluation across all 295,296 candidate pairs against empirical binary labels ($\ge 1$ dimer read).

| Model / Heuristic | Scoring Type | ROC-AUC | PR-AUC (AP) | Precision@10 | Precision@50 | Precision@100 | Brier Score |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **SADDLE Badness** | Biophysical Heuristic | 0.7166 | **0.0597** | **40.0%** | **34.0%** | **24.0%** | N/A |
| **SantaLucia ($-\Delta G$)** | Thermodynamic ($\Delta G_{37}^\circ$) | 0.7105 | 0.0125 | 20.0% | 10.0% | 6.0% | N/A |
| **Logistic Regression (Augmented)** | Supervised Probabilistic | 0.7317 | 0.0205 | 10.0% | 4.0% | 8.0% | 0.1986 |
| **XGBoost (Augmented)** | Supervised Probabilistic | **0.7564** | 0.0371 | 20.0% | 20.0% | 16.0% | 0.2143 |
| *Constant Baseline (Deterministic)* | *Deterministic Tie-Break* | *0.5000* | *0.0040* | *0.0%* | *0.0%* | *0.0%* | *0.0040* |
| *Random Choice Baseline (Expected)* | *Theoretical Expectation* | *0.5000* | *0.0040* | *0.40%* | *0.40%* | *0.40%* | *0.0040* |

: Supplementary benchmark of sequence-based heuristics and supervised models across 295,296 candidate primer pairs from the 384-plex validation panel. Supervised models were trained on 96-plex panels using thermodynamic, sequence-alignment, and approximate SADDLE-style Badness features. Binary labels represent observed dimer formation ($\ge 1$ NGS read). Precision@$K$ evaluated using deterministic secondary key sorting ($\text{rank\_key}(i) = (-\text{score}_i, \text{pair\_index}_i)$). Under this deterministic secondary sort, the Constant Baseline has zero observed dimer pairs among the first 100 pairs (yielding 0.0% at P@10, P@50, and P@100); the theoretical expected precision for an unranked random choice across the entire panel is the empirical base rate ($1,172 / 295,296 = 0.3969\% \approx 0.40\%$).

### Findings from Supplementary Benchmarking

1. **Heuristic vs. Thermodynamics:** SADDLE Badness achieves an Average Precision of 0.0597, representing a **4.8-fold improvement over nearest-neighbor free energy alone** ($-\Delta G_{37}^\circ$, AP = 0.0125) and a **4.0-fold improvement in Precision@100** (24.0% vs. 6.0%). This supports the hypothesis of Xie et al. that terminal 3' alignment and heavy GC penalties capture enzymatic extension risk more effectively than global hybridization stability.
2. **Top-Rank Discrimination vs. Supervised Models:** SADDLE Badness outperforms both Logistic Regression (AP = 0.0205, P@100 = 8.0%) and XGBoost (AP = 0.0371, P@100 = 16.0%) in top-rank precision, demonstrating that simple biophysical heuristics remain competitive at the extreme head of the distribution even when compared against supervised models incorporating Badness alongside thermodynamic parameters.
3. **Severe Probability Miscalibration of Supervised Models:** Independent recalculation of Brier scores reveals that both Logistic Regression (0.1986) and XGBoost (0.2143) are severely miscalibrated relative to the uninformative baseline ($p(1-p) = 0.003953 \approx 0.0040$). Because the supervised models were trained to maximize ranking separation under class imbalance, their mean predicted probabilities (~0.38) greatly exceed the true empirical prevalence (0.00397), yielding large quadratic calibration penalties.

# Methodological Disclosures

1. **Dual Evaluation: Primer Pairing vs. Genuine Gapped Sequence Alignment:** Rather than relying solely on primer-pairing heuristics or unaligned string comparisons, our replication evaluated both classification paradigms across all 3,942,385 reads:
   - *Method 1 (Primer Pairing):* Assigns on-target status based on cognate primer IDs ($P_2 = P_1 + 384$) and insert length threshold $\ge L_1 + L_2 - 2$ bp, yielding 43.86% on-target, 55.27% nonspecific, and 0.87% primer-dimers. This reproduces published Figure 3f (~43% / ~56% / ~1%) within a 1.0% absolute tolerance (difference = +0.86% on-target).
   - *Method 2 (Genuine Gapped Sequence Alignment via edlib):* Performs Myers' exact bit-vector semi-global sequence alignment against the assigned cognate reference amplicon from Supplementary Data MOESM5 [@xie2022designing]. At $\ge 75\%$ identity, reads partition into 43.61% on-target, 55.52% nonspecific, and 0.87% primer-dimers (difference = +0.61% on-target). At $\ge 80\%$ identity, reads partition into 42.76% on-target, 56.37% nonspecific, and 0.87% primer-dimers (difference = -0.24% on-target). Both aligner cutoffs match published Figure 3f within scientific tolerance.
   - *Positional vs Gapped Alignment Insight:* The ungapped comparison yields 39.49%, whereas gapped alignment yields 43.61% ($\ge 75\%$ edit-distance-derived similarity); a 615-read selected audit demonstrates gaps as one operational explanation, without establishing the full-library cause. While all 615 audit reads contained alignment gaps (median 8 bp), this selected sample does not establish that small 1–2 bp sequencing indels are the predominant cause of read divergence across the full 3.94-million-read library. Crucially, primer-dimer quantification (0.87% vs ~1% published) is completely invariant across all alignment methods.
   - *Alignment Scope and Concordance:* In the original SADDLE publication [@xie2022designing], on-target reads were aligned to intended amplicon sequences using Bowtie2 in MATLAB. In our lightweight, turnkey Python implementation, Method 2 executes cognate-anchored semi-global gapped alignment (`edlib`, Myers' bit-vector algorithm in prefix-matching `SHW` mode) where each cognate-primer read is aligned directly against its assigned target amplicon from MOESM5. A separate dual-strand all-target audit across an initial deterministic sample of 5,000 classified reads (768 targets per read) confirmed 0 cross-target best matches among cognate reads (0/2,203 best-scoring alignments $\ge 75\%$ were non-cognate; secondary sub-threshold alignments were not logged) and a minor 2.52% target capture among non-cognate reads, establishing 97.48% operational binary concordance across an exhaustive 2 $\times$ 2 confusion matrix ($[2,193 + 2,681] / 5,000$). Crucially, this concordance metric measures agreement specifically within this 5,000-read sample; it is not derived from a full-dataset confusion matrix, should not be extrapolated as proof of identical classification across all 3.94 million reads, and does not establish mathematical equivalence to the original Bowtie2 amplicon mapping workflow.
2. **Preprocessing Denominators:** The original publication [@xie2022designing] reported rounded percentage summaries in Figure 3f (~43% on-target, ~56% nonspecific, ~1% primer-dimers) without disclosing the underlying integer read counts or their internal filtering denominator. While our pipeline's derived post-filter classified denominator ($N = 2,860,664$) yields proportions in close numerical proximity (43.61% / 55.52% / 0.87%), the original authors' internal denominator remains unverified. We therefore report both our derived post-filter denominator ($N = 2,860,664$) and the total raw sequencer output ($N = 3,942,385$, across which on-target reads represent 31.64%), presenting agreement strictly as approximate operational concordance under our documented filtering criteria rather than an established identical denominator.
3. **Qualified Partial Replication Scope Limitation:** This study constitutes an independent **qualified partial replication** focused strictly on the empirical validation read partitions and heuristic scoring. SADDLE's simulated annealing primer pool design software was withheld under patent protection and is not evaluated. In addition, the read classification is reimplemented in Python using cognate-anchored semi-global gapped alignment against target amplicons rather than Bowtie2 amplicon mapping in MATLAB.
4. **Distributional Terminology:** We characterize the distribution of primer-dimer read counts as heavy-tailed and highly concentrated; no formal power-law distribution fit (e.g., Clauset-Shalizi-Newman methodology) was conducted.
5. **Data and Code Licensing:** Raw sequencing data from Figshare (DOI: 10.6084/m9.figshare.16944154.v3) are licensed under CC BY 4.0. The original analysis scripts were made publicly available through GitHub and Figshare without an explicit license header. The replication codebase is distributed under the MIT License.

# Conclusion

We have performed an independent **qualified partial computational replication** of the empirical NGS validation results from SADDLE [@xie2022designing]. Using an open-source Python workflow, we confirmed that classified reads from the SADDLE-designed 384-plex validation panel (Figure 3f) replicate published values (~43% on-target, ~56% nonspecific, ~1% primer-dimers). On our pipeline's derived post-filter classified denominator ($N = 2,860,664$), reads partition under primer-pairing classification (Method 1) into 43.86% on-target, 55.27% nonspecific, and 0.87% primer-dimers. Under genuine gapped sequence alignment against assigned target amplicons (Method 2 via `edlib`), reads partition into 43.61% on-target ($\ge 75\%$ identity) or 42.76% on-target ($\ge 80\%$ identity), 55.52% to 56.37% nonspecific, and 0.87% primer-dimers (24,878 reads), approximately reproducing the published ~43% / ~56% / ~1% figures within 0.61% to 0.86% absolute difference. Evaluated across all 3,942,385 raw sequencing reads, on-target amplicons represent 31.64% of sequencer output. Under all criteria and denominators, the observed primer-dimer fraction remains low (0.87% of classified reads, 0.63% of raw reads), which is consistent with the published experimental finding that the SADDLE-designed panel exhibited suppressed dimer formation. We also corroborated the severe heavy-tailed concentration of dimer reads, with a single primer pair generating 67.57% of all dimer reads in the verified pair-level distribution (16,805 / 24,872 reads in `data/saddle_384_dimer_pairs.tsv`). Finally, our supplementary exploratory benchmark illustrated the predictive utility of the SADDLE Badness heuristic over standard thermodynamics while documenting the calibration characteristics of supervised models.

All data, scripts, and a standalone verification suite are available in the public replication repository.

# Code and Data Availability

- **GitHub Repository:** [https://github.com/Rashidmstar12/saddle-replication](https://github.com/Rashidmstar12/saddle-replication)
- **Zenodo Deposit:** [@kadayil2026saddle, https://doi.org/10.5281/zenodo.24100087]
- **Primary Data Source:** Figshare (DOI: [10.6084/m9.figshare.16944154.v3](https://doi.org/10.6084/m9.figshare.16944154.v3), CC BY 4.0)

Reviewers can verify all archived primary and supplementary results, benchmark tables, and member file checksums in under two seconds using the standalone verifier (full raw FASTQ processing requires running `reproduce_from_fastq.py --full` with external raw sequencing inputs):

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
