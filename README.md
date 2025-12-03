# cfRNA-Microbiome

Code and intermediate data for the manuscript:

This repository contains the scripts and processed tables used to reproduce the main figures and summary statistics in the manuscript. It focuses on **host-filtered blood RNA sequencing data** and downstream **Kraken2 / MetaPhlAn4**-based taxonomic profiling.

---

## Contents

- `Script/`  
  Shell and R scripts for:
  - host read removal and BAM/FASTQ processing  
  - non-host read extraction  
  - Kraken2 / MetaPhlAn4 classification  
  - aggregation of per-sample taxonomic profiles  
  - generation of summary tables and figures (including those used in the paper)

  See the header comments of each script for input / output paths and usage.

- `intermediate_data/`  
  Processed data tables used directly by the R visualization scripts, for example:
  - per-sample read counts and host/non-host fractions  
  - per-cohort taxonomic abundance matrices  
  - summary statistics used in Figures and Tables S1–S2  

- `LICENSE`  
  MIT license for this repository.

---

## Data availability

This repository does **not** redistribute raw sequencing data.

- CAD whole-blood RNA-seq: GEO **GSE58150**  
- Tuberculosis plasma cfRNA: GEO **GSE255073**

Please download raw data from GEO and perform initial alignment/host-filtering in your own environment. The scripts here assume that host-filtered BAM/FASTQ files and taxonomic classification outputs are already available.

Processed data tables and analysis code in this repository correspond to the materials described in the manuscript and its Supplementary Tables.

---

## Reproducing analyses and figures

1. **Prepare inputs**

   - Download raw FASTQ files from GEO (GSE58150, GSE255073).  
   - Perform host alignment and filtering (e.g., Bowtie2 + SAMtools) to obtain non-host read sets.  
   - Run Kraken2 and/or MetaPhlAn4 with the reference databases used in the manuscript.  

2. **Generate intermediate tables**

   - Use the shell and R scripts in `Script/` to parse classifier outputs, aggregate read counts, and build per-sample abundance matrices.  
   - Output files should be written into `intermediate_data/` as indicated in each script.

3. **Visualize results**

   - Run the R visualization scripts (e.g. `05_results_visualization.R`) to reproduce summary plots and tables for the CAD and TB cohorts.

Because paths and module systems are environment-specific, you will likely need to adjust file paths, module load commands, and number of threads in the shell scripts.

---

## Software requirements

The analysis was developed and run on a Linux HPC environment. In general you will need:

- A recent version of **R** (≥ 4.x) with standard tidyverse/plotting packages  
- Command-line tools for alignment and metagenomics (e.g. Bowtie2, SAMtools, Kraken2, MetaPhlAn4)  

See individual script headers for the exact package and tool requirements.

---

## Citation

If you use this code or processed data in your own work, please cite the corresponding manuscript:

> Kui Chen *et al.* Host-filtered Blood Nucleic Acids for Pathogen Detection: Shared Background, Sparse Signal, and Methodological Limits. *Pathogens* (2025), in review.

---

## Contact

For questions about the code or analyses, please contact:

Zhaoxia Wang 
wangzhaoxia1004@126.com

Issues and pull requests are welcome.
