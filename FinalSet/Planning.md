

> **Project:** Benchmarking detection of DNA copy number variation using high-density methylation microarrays 
> **Lab:** Przybyl Lab · RI-MUHC / McGill University 
> **Author:** Eliza Kishan 
> **Data source:** `FinalSet`

---

## Central Thesis

> The accuracy of methylation-based CNV calling tools is not determined by software choice or parameter tuning, it is determined by the **genomic complexity of the tumour**, as quantified by the Genomic Index (GI).
> 
> In low-complexity benign tumours (LM), all three tools perform adequately. In high-complexity malignant tumours (LMS), all three tools fail **proportionally** to the degree of genomic rearrangement. 

### At-a-glance stats

|Metric|Value|
|---|---|
|LMS mean accuracy (all tools)|53–58%|
|LM mean accuracy (all tools)|92–96%|
|GI–accuracy Pearson r|−0.66|
|GI seen by tools vs SNP|13–22% of true GI|

---

## 1. Introduction

### Structure

- LMS and LM are smooth-muscle tumours with opposite genomic profiles. LMS is malignant and has high genomic complexity; LM is benign and has low genomic complexity.
- SNP arrays are the gold standard for CNV detection but require fresh tissue and specialised infrastructure. Methylation arrays (EPIC) are increasingly used for CNV calling because the data is already being generated for epigenetic profiling. No benchmarking study has systematically compared the two platforms across the full benign-to-malignant spectrum.
- **New Contribution:**
    1. Accuracy benchmarking across three tools and four bin sizes
    2. A mechanistic explanation via the Genomic Index
    3. Gene-level concordance for clinically relevant loci
    4. A null result on bin size that rules out a parameter explanation for observed failures
- The accuracy numbers reported are concordant with an [existing paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7028453/). The novel contribution is not the accuracy numbers themselves but that the scope of use might be limited by the genomic complexity.

### Key citations needed

- CBS (circular binary segmentation) algorithm 
- Sesame, MethylMaster, and Conumee tool papers
- Concordant accuracy paper — introduce here, revisit in Discussion
- GI as a measure of genomic complexity — define formally here or in Methods

---

## 2. Methods

### What to include

- **Dataset:** 13 LMS cases and 9 LM cases. All processed via the EPIC methylation array. SNP data from OncoScan CNV arrays (gold standard). 
- **Tools evaluated:** Sesame, MethylMaster, Conumee, each at default + four bin sizes (10 kb, 50 kb default, 100 kb, 1 Mb).
- **Preprocessing:** Raw data converted to standardised GenomicRanges dataframes.
- **Genomic Index:** It penalises both the number and amplitude of copy number events. Formula below.
$$GI = \frac{(\sum \# \text{CNVs})^2}{\text{\# of chromosomes showing somatic changes}}$$
- **Accuracy metric 1 — whole-dataset:** Overlapping predicted and truth GenomicRanges, matching basepairs divided by total genome assessed. Mean taken across all cases.
- **Accuracy metric 2 — median:** Same calculation, but median basepair accuracy per case. 
- **Gene concordance:** Seven genes assessed: CDKN2A, CCNE1, MYC, MYOCD, PTEN, RB1, TP53. CNV status from each methylation tool vs SNP ground truth.
- **Genome modified:** Percentage of genome covered by called CNV segments per case and tool. Assesses detection sensitivity independent of call accuracy.

---

## 3. Results

### 3.1 Accuracy splits cleanly by tumour type
#### Numbers to cite

|Tool|LMS mean accuracy|LM mean accuracy|
|---|---|---|
|Conumee|58.4%|91.6%|
|MethylMaster|54.8%|**96.4%**|
|Sesame|52.8%|95.9%|

- **LM median accuracy >99.8%** for all tools. LM failures are concentrated in a small number of outlier cases (primarily Case 9312, LM outlier), not a systematic failure.
- **Tool ranking flips between cohorts.** Conumee leads in LMS but is worst in LM. MethylMaster leads in LM. 

#### Figures

- **Fig. 5 (poster)** — Bin-wise accuracy bar chart, LMS. Point out that cases 9338, 9328, 9350 cluster at the bottom regardless of tool or bin.
- **Fig. 6 (poster)** — Bin-wise accuracy bar chart, LM. Highlight consistency, with 9312 as the visible exception.

---

### 3.2 Bin size does not rescue accuracy in LMS

**A null result. **

|Bin size|LMS mean accuracy|Note|
|---|---|---|
|10 kb|54.2%|Baseline|
|50 kb|55.3%|No improvement|
|100 kb|55.2%|No improvement|
|1 Mb|56.7%|Marginal — not meaningful|

Accuracy is flat across two orders of magnitude. If CBS failure were a resolution problem, you would expect a clear recovery at 1 Mb. The 2.5 pp increase at 1 Mb is not meaningful. 

**This means that there is no explanation on the basis of the bin parameter.**

#### Figures

- **Fig. 5 & 6 (poster)** — Highlight in text that bars are nearly identical across bins within each case.

---

### 3.3 Genomic Index explains the accuracy gradient

**Argument built in 3 steps: **

1. GI separates cases by complexity
2. Accuracy decreases with GI
3. Methylation tools systematically underestimate GI.

#### Case-level LMS table (source: `GenomicIndexLMS.csv` + `accuraciesLMS.csv`)

| Case     | SNP GI    | SNP genome % | Mean accuracy | Flag                                  |
| -------- | --------- | ------------ | ------------- | ------------------------------------- |
| 9350     | 3,463     | 40.9%        | 31.6%         | Highest GI in cohort, re-verified     |
| 9328     | 3,388     | 22.3%        | 28.1%         | Extreme GI, re-verified               |
| 9338     | 3,002     | 16.8%        | 24.2%         | Worst overall, re-verified            |
| 9357     | 917       | 37.7%        | 54.7%         | Average                               |
| 9356     | 666       | 41.5%        | 56.0%         | Average                               |
| 9353     | 550       | 47.3%        | 48.1%         | Average                               |
| 9355     | 511       | 42.8%        | 53.3%         | Average                               |
| 9354     | 393       | 45.7%        | 56.0%         | Average                               |
| 9327     | 464       | 42.1%        | 52.4%         | Average                               |
| 9358     | 868       | 44.3%        | 45.0%         | Average                               |
| **9337** | **2,200** | **25.1%**    | **77.7%**     | **ANOMALY**                           |
| 9202     | 60        | 16.7%        | 94.3%         | LM-equivalent GI, re-verified by hand |
| 9203     | 56        | 21.4%        | 98.1%         | LM-equivalent GI, re-verified by hand |

**Pearson r = −0.659** (SNP GI vs mean accuracy across LMS cases) 
**Pearson r = −0.652** (SNP GI vs mean accuracy across LM cases)
***This shows decent correlation.***
#### GI compression ratios (source: `GenomicIndexLMS.csv` + `genomeChangedLMS.csv`)

| Tool         | GI seen as % of SNP GI | Genome territory detected  | Failure mode                    |
| ------------ | ---------------------- | -------------------------- | ------------------------------- |
| MethylMaster | 21.6%                  | **96.4% of SNP territory** | Finds regions; misses CNV depth |
| Conumee      | 16.2%                  | 76.1% of SNP territory     | Partial territory + depth loss  |
| Sesame       | 12.8%                  | **49.2% of SNP territory** | Misses half the affected genome |

> [!important] Two mechanistically distinct failure modes: **MethylMaster** correctly identifies 96% of the genome alteration but compresses copy number amplitude. **Sesame** fails at the location level, detecting only 49% of affected genome territory, and also compresses amplitude. This can be verified with a manual observation of the graphs generated.

#### Case 9337, anomalous

> [!note] **SNP GI = 2,200**, but achieves **77.7% mean accuracy**.
> 
> Methylation tools detect _more_ genome territory than SNP does in this case (MethylMaster ratio = 1.24×, Conumee = 1.26×). Manual observation showed that the high GI was actually driven by a sort of extreme segmentation in the SNP data, which means that the value was extremely inflated, **HOWEVER**, notably, they were in widely segmented, so by just measuring deletions and amplifications, it gives a higher accuracy.
> 
> **Implication:** GI is a necessary but not sufficient predictor of tool performance. 

#### Cases 9202 and 9203

SNP GIs of 60 and 56, close to the LM cohort (median SNP GI ~8). Accuracies of 94.3% and 98.1% track accordingly.

**Including them in LMS summary stats artificially elevates the mean.** 
*Options to mitigate statistical inaccuracy:*
- Report LMS accuracy with and without these two cases
- Flag them as a biologically distinct sub-group of low-complexity LMS
#### Figures
- **Fig. 7 (poster)** — LMS and LM metrics heatmap. 

---

### 3.4 Gene concordance reveals selective blind spots

#### Gene concordance table

| Gene      | Concordance | SNP ground truth    | Interpretation              |
| --------- | ----------- | ------------------- | --------------------------- |
| CDKN2A    | **87%** ✅   | Deleted 12/13 LMS   | Large canonical deletion    |
| TP53      | **85%** ✅   | Amplified 9/13 LMS  | Strong amplification signal |
| PTEN      | 74%         | Mixed               | Acceptable                  |
| RB1       | 62%         | Mixed               | Borderline                  |
| MYC       | 41% ❌       | Focal amplification | Near-random                 |
| CCNE1     | 41% ❌       | Focal amplification | Near-random                 |
| **MYOCD** | **36%** ❌   | Focal               | **Most significant miss**   |

- **Deletion vs amplification:** CDKN2A (deletion, 87%) and TP53 (amplification, 85%) are captured with similar success. This means that the distinction between the difficulty of calling amplification v. deletion might not necessarily be significantly different.

#### Figures

- **Fig. 4 (poster)** — Gene concordance heatmap. 

---

## 4. Discussion

### Outline

**¶1 — Central finding:** Methylation-based CNV calling in LMS performs near chance (~53–58%), while LM performance is clinically adequate (~92–96%). 

**¶2 — Connect to the concordant prior paper:** Your accuracy numbers reproduce those reported by [Smith et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7028453/). 

**¶3 — The CBS failure mechanism:** CBS segments by identifying statistically significant transitions in log2 copy number ratios. In a genome where 25–47% of chromosomal territory is altered (mid-range LMS cases), the noise floor is elevated globally. The null result on bin size (flat from 10 kb to 1 Mb) confirms this is not a resolution artefact.

**¶4 — Two failure modes of tools:** MethylMaster and Sesame represent mechanistically distinct failure modes:

- MethylMaster: high sensitivity to the genome modified (96% territory), compresses amplitude
- Sesame: misses 51% of genome modified, also compresses amplitude

Future improvement efforts should target different aspects of each tool's pipeline.

**¶5 — The 9337 anomaly and CNV architecture:** Case 9337 demonstrates that aggregate GI is not the only predictor. Despite SNP GI of 2,200, it achieves 77.7% accuracy. Genome-changed analysis shows methylation tools detect more chromosomal territory than SNP in this case (ratios >1.2×). Manual verification shows that the modification in the genome of case 9337 is focal, but just that the segmentation of the SNP data seems to be very fine, leading to an inflated GI. This however does show that the structure of the CNVs is also important to accuracy as much as the algorithm is.

**¶6 — Gene discovery implications:** The methylation platform reliably detects large canonical events (CDKN2A, TP53) but fails on focal events (MYOCD, CCNE1, MYC). MYOCD concordance of 36% is particularly concerning given its biological relevance in LMS. 

**¶7 — Recommendations**

| Cohort          | Recommendation                                                                                                                         |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| LM (benign)     | All three tools clinically adequate. MethylMaster/Sesame preferred. Can be used without SNP validation for most research applications. |
| LMS (malignant) | Should not replace SNP arrays for individual case-level  decisions. Appropriate for large-scale population-level analyses.             |
| Tool selection  | No universal winner. Conumee edges out in LMS; MethylMaster/Sesame lead in LM.                                                         |
| GI threshold    | Cases with high GI might be high-risk for low methylation-platform accuracy.                                                           |

**¶8 — Limitations**

- Small cohort: 13 LMS and 9 LM cases. Case-level conclusions (9337, 9202, 9203) are hypothesis-generating.
- Single platform: all methylation data from EPIC array. Used for normalization and to prevent bias from different platforms, performance may differ on 450k (might be contraindicated by the Smith paper) or next-generation arrays.
- No CBS alternatives tested: SLM, GISTIC, or HMM-based segmentation are future directions.

---

## 5. Conclusion
1. **Finding:** Methylation-based CNV calling achieves high accuracy in benign LM tumours (92–96%) but performs near chance in malignant LMS (53–58%).
2. **Mechanism:** This gap is explained by the Genomic Index, accuracy is inversely correlated with genomic complexity (r = −0.66), and is not sensitive to bin size changes.
3. **Gene level:** Large canonical events (CDKN2A, TP53) are reliably detected, while focal events (MYOCD, CCNE1, MYC) are not.
4. **Implication:** These findings define the conditions under which methylation-based CNV calling is accurate and posit that GI is a major condition in its performance.

---

## 6. Figure List

### Figures already produced and in poster

| Figure | Description                                                  | Where to use       |
| ------ | ------------------------------------------------------------ | ------------------ |
| Fig. 1 | Genomic Index method diagram                                 | Methods            |
| Fig. 2 | Gene concordance method diagram                              | Methods            |
| Fig. 3 | Whole-dataset accuracy method diagram                        | Methods            |
| Fig. 4 | Gene concordance heatmap — 7 genes × 12 method/bin combos    | Results §3.4       |
| Fig. 5 | Bin-wise accuracy bar chart — LMS cases                      | Results §3.1, §3.2 |
| Fig. 6 | Bin-wise accuracy bar chart — LM cases                       | Results §3.1       |
| Fig. 7 | LMS and LM metrics heatmap — clustered by GI and performance | Results §3.3       |
### Unused Figures (Not in Poster)

| Figure                                         | Data needed                                    | Key message                           |
| ---------------------------------------------- | ---------------------------------------------- | ------------------------------------- |
| Per-case segment plots comparing tools and SNP | `GenomicIndexLMS.csv`                          | compression of amplitude by all tools |
| Side-by-side: LMS vs LM mean accuracy by tool  | `accuraciesLMS.csv` + `accuraciesLM.csv`       |                                       |
| Genome Changed Bar Graph                       | `genomeChangedLMS.csv` + `genomeChangedLM.csv` | Null result, flat trajectories       |

### Potential figures to create and add

| Figure                                                  | Data needed                                    | Key message                                          |
| ------------------------------------------------------- | ---------------------------------------------- | ---------------------------------------------------- |
| Scatter: SNP GI vs mean accuracy (LMS cases, annotated) | `GenomicIndexLMS.csv` + `accuraciesLMS.csv`    | r = −0.659; annotate 9337, 9202/9203, 9338/9328/9350 |
| Genome-changed violin: LMS vs LM per tool               | `genomeChangedLMS.csv` + `genomeChangedLM.csv` | MethylMaster 96% territory vs Sesame 49%             |
| Accuracy line chart: bins 10 kb→1 Mb per tool (LMS)     | `accuraciesLMS.csv`                            | Null result, flat trajectories                      |

---

## 7. Statistics Reference

### Accuracy

|Metric|Value|Source|
|---|---|---|
|LMS mean accuracy — Conumee|58.4%|`accuraciesLMS.csv`|
|LMS mean accuracy — MethylMaster|54.8%|`accuraciesLMS.csv`|
|LMS mean accuracy — Sesame|52.8%|`accuraciesLMS.csv`|
|LM mean accuracy — MethylMaster|96.4%|`accuraciesLM.csv`|
|LM mean accuracy — Sesame|95.9%|`accuraciesLM.csv`|
|LM mean accuracy — Conumee|91.6%|`accuraciesLM.csv`|
|LM median accuracy (all tools)|>99.8%|`accuraciesMedianLM.csv`|
|LMS accuracy range by bin (all tools)|54.2–56.7%|`accuraciesLMS.csv`|

### Genomic Index

|Metric|Value|Source|
|---|---|---|
|Mean SNP GI — LMS cohort|~1,272 (default bin)|`GenomicIndexLMS.csv`|
|Mean SNP GI — LM cohort|~14 (default bin)|`LMGenomicIndex.csv`|
|Highest SNP GI — case 9350|3,463|`GenomicIndexLMS.csv`|
|GI compression — MethylMaster|21.6% of SNP GI|`GenomicIndexLMS.csv`|
|GI compression — Conumee|16.2% of SNP GI|`GenomicIndexLMS.csv`|
|GI compression — Sesame|12.8% of SNP GI|`GenomicIndexLMS.csv`|
|GI–accuracy Pearson r (LMS)|−0.659|Derived|
|GI–accuracy Pearson r (LM)|−0.652|Derived|

### Genome Changed

|Metric|Value|Source|
|---|---|---|
|Mean genome changed — SNP (LMS)|~34% across cases|`genomeChangedLMS.csv`|
|MethylMaster detection rate|96.4% of SNP territory|`genomeChangedLMS.csv`|
|Conumee detection rate|76.1% of SNP territory|`genomeChangedLMS.csv`|
|Sesame detection rate|49.2% of SNP territory|`genomeChangedLMS.csv`|
|9337 MethylMaster ratio|1.24× (detects more than SNP)|`genomeChangedLMS.csv`|
|9337 Conumee ratio|1.26× (detects more than SNP)|`genomeChangedLMS.csv`|

### Gene Concordance

| Gene      | Concordance | SNP status (LMS)    | Clinical relevance    |
| --------- | ----------- | ------------------- | --------------------- |
| CDKN2A    | **87%**     | Deleted 12/13       | Tumour suppressor     |
| TP53      | **85%**     | Amplified 9/13      | Tumour suppressor     |
| PTEN      | 74%         | Mixed               | Tumour suppressor     |
| RB1       | 62%         | Mixed               | Cell cycle regulation |
| MYC       | 41%         | Focal amplification | Oncogene              |
| CCNE1     | 41%         | Focal amplification | Cell cycle            |
| **MYOCD** | **36%**     | Focal               | **LMS-specific**      |

---
All data can be found in FinalSet for reverification.
