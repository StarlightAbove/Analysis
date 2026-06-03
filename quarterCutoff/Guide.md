---
title: "Guide"
output: html_document
---

# quarterCutoff — Analysis Directory

This directory contains the outputs and source code for the quarter-cutoff CNV analysis, comparing leiomyoma (LM) and leiomyosarcoma (LMS) samples against matched normals.

------------------------------------------------------------------------

## Directory Structure

### `analyticFile.R`

The master analysis script. All code is organized under named headings, covering each analytical module described below.

------------------------------------------------------------------------

### `case_graphs/`

Per-sample genome-wide copy number plots in PDF format. Each sample is represented at four bin resolutions: **10,000**, **50,000**, **100,000**, and **1,000,000** bp.

Samples are organized into three subfolders:

-   **`LM/`** — Leiomyoma cases (IDs: 9312, 9314, 9316, 9320, 9322, 10067, 10069, 10090, 10092). Plots show genome-wide copy number profiles.
-   **`LMS/`** — Leiomyosarcoma cases (IDs: 9202, 9203, 9328, 9337, 9338, 9350, 9353–9358). Plots additionally show log₂ ratios for key genes of interest.
-   **`Normals/`** — Matched normal controls (Illumina array IDs). Used as the reference baseline for copy number calculations.

------------------------------------------------------------------------

### `Accuracy/`

Genome-wide accuracy metrics comparing the CNV calls against the ground truth.

| File              | Description                     |
|-------------------|---------------------------------|
| `accuracyLM.csv`  | Accuracy values for LM samples  |
| `accuracyLM.png`  | Accuracy plot for LM samples    |
| `accuracyLMS.csv` | Accuracy values for LMS samples |
| `accuracyLMS.png` | Accuracy plot for LMS samples   |

------------------------------------------------------------------------

### `Gene_concordance/`

Concordance between CNV calls from the tested technology and the ground truth, assessed per gene.

| File | Description |
|----|----|
| `CCNE1.csv`, `CDKN2A.csv`, `MYC.csv`, `MYOCD.csv`, `PTEN.csv`, `RB1.csv`, `TP53.csv` | Per-gene concordance values |
| `geneConcordance.csv` | Aggregated concordance table across all genes |
| `geneConcordance.png` | Summary plot for an at-a-glance view of concordance across genes |

------------------------------------------------------------------------

### `Genomic_Index/`

Genomic index (GI) calculations, used as a measure of genomic instability.

| File                  | Description                |
|-----------------------|----------------------------|
| `genomicIndexLM.csv`  | GI values for LM samples   |
| `genomicIndexLMS.csv` | GI values for LMS samples  |
| `LMGI.png`            | GI comparison plot for LM  |
| `LMSPlot.png`         | GI comparison plot for LMS |

------------------------------------------------------------------------

### `Genome_modified/`

Proportion of the genome showing copy number alterations.

| File                      | Description                               |
|---------------------------|-------------------------------------------|
| `genome_modifiedLM.csv`   | Genome-modified fractions for LM samples  |
| `genome_modified_LMS.csv` | Genome-modified fractions for LMS samples |
| `GenomeModifiedLM.png`    | Comparison plot for LM                    |
| `genomeModifedLMS.png`    | Comparison plot for LMS                   |

------------------------------------------------------------------------

### `normals.csv`

Summary table for the normal control samples used in this analysis.
