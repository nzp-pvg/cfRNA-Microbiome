# 🧬 cfRNA-Microbiome Reproducibility Hub

> **Host-filtered Blood Nucleic Acids for Pathogen Detection: Shared Background, Sparse Signal, and Methodological Limits**  
> Re-analysis of host-filtered whole-blood RNA-seq and plasma cfRNA data from CAD and tuberculosis cohorts.

[![Analysis](https://img.shields.io/badge/analysis-reproducible-0F766E?style=for-the-badge)](#-reproducibility-workflow)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0-276DC3?style=for-the-badge&logo=r)](#-software-requirements)
[![Shell](https://img.shields.io/badge/shell-Bash-4EAA25?style=for-the-badge&logo=gnu-bash&logoColor=white)](#-software-requirements)
[![Kraken2](https://img.shields.io/badge/taxonomy-Kraken2-8B5CF6?style=for-the-badge)](#-reproducibility-workflow)
[![MetaPhlAn4](https://img.shields.io/badge/taxonomy-MetaPhlAn4-E05A33?style=for-the-badge)](#-reproducibility-workflow)

## 🧭 Project Overview

This repository contains the analysis scripts, intermediate taxonomic profiles, summary tables, manuscript materials, and submission records used to examine whether host-filtered blood nucleic-acid sequencing contains a reproducible microbial signal. The workflow is designed to separate **technical background**, **classification sensitivity**, and **disease-associated signal** rather than treating every computationally classified read as evidence of infection.

The analysis combines:

- 🧹 SRA download, FASTQ conversion, compression, and `fastp` quality control
- 🧬 Host-filtering sensitivity analysis using **strict** and **relaxed** settings
- 🦠 Kraken2 report parsing and classified/unclassified read summaries
- 🌿 MetaPhlAn4 profile aggregation at species level
- 📊 Cross-cohort comparison of CAD and tuberculosis samples
- 🧪 Exploratory assessment of *Mycobacterium tuberculosis* signal and Kraken2/MetaPhlAn4 concordance
- 📝 Generation of figure-source tables, supplementary summaries, and manuscript-ready materials

## 🎯 Scientific Focus

The central analytical question is:

> **After host-read removal, how much of the remaining blood RNA-seq signal is reproducible microbial background, and how much can be interpreted as disease-associated evidence?**

The repository therefore emphasizes three safeguards:

1. **Filtering sensitivity:** compare strict and relaxed host-filtering modes.
2. **Method concordance:** compare Kraken2 and MetaPhlAn4 where both profiles are available.
3. **Interpretation discipline:** computational taxonomic assignments are treated as signals requiring validation, not as direct proof of viable pathogens or causality.

## 🧬 Cohorts and Data Availability

Raw sequencing data are **not redistributed** in this repository. Users should obtain the source data from the relevant public repositories and comply with their access and reuse conditions.

| Cohort | Accession | Material / analysis role |
|---|---|---|
| CAD | [GSE58150](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58150) | Whole-blood RNA-seq; cardiovascular disease and control comparison |
| Tuberculosis | [GSE255073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255073) | Plasma cfRNA; tuberculosis-positive and tuberculosis-negative comparison |

The repository includes processed Kraken2 reports and MetaPhlAn4 profiles used for the reported analyses. Raw FASTQ files, host-reference indices, Kraken2 databases, MetaPhlAn4 databases, and cluster-specific temporary files must be prepared separately.

## 🔬 Reproducibility Workflow

```mermaid
flowchart LR
    A[Public sequencing accessions] --> B[SRA Toolkit]
    B --> C[FASTQ conversion]
    C --> D[fastp quality control]
    D --> E[Host filtering]
    E --> F{Filtering mode}
    F --> G[Strict non-host reads]
    F --> H[Relaxed non-host reads]
    G --> I[Kraken2 reports]
    H --> I
    G --> J[MetaPhlAn4 profiles]
    H --> J
    I --> K[Read-level QC summaries]
    J --> L[Species-level abundance summaries]
    K --> M[Figures and supplementary tables]
    L --> M
```

### 1. 📥 Input preparation

`int/Script/01_download_raw_fastq.sh` downloads SRA records from an accession list, converts them to paired-end FASTQ, compresses the files, and runs `fastp`. Before execution, update `PROJECT_ROOT`, `ACC_LIST`, and the output directories in the script.

### 2. 🧹 Host filtering and taxonomic classification

The downstream scripts expect host-filtered outputs and classifier results arranged by cohort, filtering mode, and biological group. The current project stores these results under `int/res/`.

Expected group labels are:

- `GSE58150`: `ctl` and `cvd`
- `GSE255073`: `neg` and `pos`

Expected filtering modes are `strict` and `relaxed`.

### 3. 🦠 Kraken2 summaries

The Kraken2 workflow reads `.kreport` files and derives:

- total classified and unclassified reads
- classified and unclassified percentages
- group-by-mode summaries
- per-sample strict-minus-relaxed differences
- Welch *t*-tests and permutation-based comparisons
- strict-mode summaries and Prism-ready figure-source tables

The main analysis code is located at `int/res/ms_5.R`. Because this file contains legacy analysis blocks from the manuscript development process, inspect and update the input/output paths before running it as a clean end-to-end pipeline.

### 4. 🌿 MetaPhlAn4 summaries

`int/Script/cvd_ms_5_res_MetaPhlAn.R` reads MetaPhlAn4 profiles, extracts species-level abundances, calculates classified fraction and richness, and creates summary plots including:

- classified percentage by filtering mode
- detected species richness
- strict-mode top-species heatmap

The script also contains optional code for exporting plot files and should be run after the MetaPhlAn4 profile directories are prepared.

## 🗂️ Repository Layout

```text
cfRNA-Microbiome/
├── Close/                         # Accepted manuscript and final publication materials
│   ├── Supplementary_Table_S1_S2_S3.zip
│   ├── pathogens-4053565_v3.docx
│   └── pathogens-4053565_v3.pdf
├── int/                            # Intermediate analysis and manuscript-development materials
│   ├── Script/                     # SRA/FASTQ and MetaPhlAn4 scripts
│   ├── fig/                        # Figure and table materials
│   ├── ref/                        # Reference materials
│   └── res/                        # Kraken2 reports, MetaPhlAn4 profiles, and summaries
├── revision/                       # Revision files, response letters, and reviewer materials
├── submission/                     # Submission package and supplementary files
└── README.md                       # Project documentation
```

## 📦 Key Intermediate Outputs

| Output type | Location / example |
|---|---|
| Kraken2 reports | `int/res/GSE58150/Kraken2_kreport/` and `int/res/GSE255073/Kraken2_kreport/` |
| MetaPhlAn4 profiles | `int/res/*/MetaPhlAn_profiles/metaphlan/` |
| MetaPhlAn4 summaries | `int/res/*/MetaPhlAn_profiles/summary_metaphlan/` |
| Manuscript figures and tables | `int/fig/Figures_Tables/` |
| Supplementary tables | `Close/Supplementary_Table_S1_S2_S3.zip` and `submission/` |
| Manuscript and submission records | `Close/`, `revision/`, and `submission/` |

## 🚀 Getting Started

### Prerequisites

Install or load the following tools in a Linux/HPC environment or an equivalent local setup:

- **R** ≥ 4.0 with `tidyverse`, `dplyr`, `readr`, `stringr`, `tidyr`, `purrr`, `writexl`, `igraph`, `tidygraph`, `ggraph`, and `scales`
- **SRA Toolkit**: `prefetch`, `fasterq-dump`
- **fastp** and `gzip`
- **Bowtie2** and **SAMtools** or the host-filtering tools used for the original analysis
- **Kraken2** with the database version used in the study
- **MetaPhlAn4** with the database version used in the study

### Suggested execution order

```bash
# 1. Prepare an accession list and edit paths in the download script
bash int/Script/01_download_raw_fastq.sh

# 2. Perform host filtering and generate Kraken2 / MetaPhlAn4 outputs
#    This step depends on the local reference databases and HPC environment.

# 3. Update ROOT / RES_ROOT paths in the R scripts, then run the summaries
Rscript int/Script/cvd_ms_5_res_MetaPhlAn.R
Rscript int/res/ms_5.R
```

⚠️ **Path note:** several legacy script blocks still contain paths from the earlier `CVD_MS_5` project layout. Change those paths to the corresponding `CVD_MS_4/int/` locations before execution. The repository preserves these scripts and outputs for traceability; a clean rerun should use explicitly versioned input and output directories.

## ♻️ Reproducibility Boundaries

- Results depend on the exact host reference, Kraken2 database, MetaPhlAn4 database, software versions, and filtering parameters.
- Strict and relaxed filtering modes are analytical conditions, not interchangeable biological measurements.
- Taxonomic classification of host-filtered blood nucleic acids is vulnerable to low biomass, index hopping, environmental background, database incompleteness, and read-level misclassification.
- Kraken2 and MetaPhlAn4 outputs should be interpreted as computational evidence and require orthogonal experimental validation.
- The archived manuscript, revision, and submission folders are retained to preserve the project history. Do not overwrite these historical materials when generating new versions.

## 📚 Citation

If you use this repository, scripts, or processed tables, please cite the associated manuscript:

**https://doi.org/10.3390/pathogens15010055**

Please also cite the original data records and the software/database resources used in your rerun.

## 📄 License

See `LICENSE` if present. If no license file is included in a public release, the code and processed materials should be treated as **all rights reserved** until an explicit license is added by the authors.

## 💬 Contact and Contributions

For questions about the analysis or manuscript materials, please contact the corresponding authors listed in the manuscript. Issues and pull requests are welcome for documentation fixes, reproducibility improvements, and clearly documented methodological corrections.

<p align="center">
  🧬 Host-filtered signal&nbsp;&nbsp;•&nbsp;&nbsp;🦠 Taxonomic profiling&nbsp;&nbsp;•&nbsp;&nbsp;📊 Transparent comparison&nbsp;&nbsp;•&nbsp;&nbsp;♻️ Reproducible analysis
</p>
