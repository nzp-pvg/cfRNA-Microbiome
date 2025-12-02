#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Script: 02_host_mapping_and_nonhost_extraction.sh
# Author: Zhaoxia Wang
# Version: v1.0
# Date: 2024-09-30
#
# Purpose:
#   1) Map plasma cfRNA FASTQ files to the human reference genome (GRCh38)
#      using Bowtie2.
#   2) Extract non-host reads under two host-filtering modes:
#        - strict mode:    paired-unmapped (both mates unmapped)
#        - relaxed mode:   any-unmapped   (at least one mate unmapped)
#
#   The script is written to support both CAD and TB cohorts, assuming
#   standard paired-end FASTQ naming: SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz
#
# Usage (example):
#   bash 02_host_mapping_and_nonhost_extraction.sh CAD
#   bash 02_host_mapping_and_nonhost_extraction.sh TB
#
# Required edits before running:
#   1) Set PROJECT_ROOT to your project directory.
#   2) Set BT2_INDEX_PREFIX to the prefix of your Bowtie2 GRCh38 index.
###############################################################################

# ------------------------- USER CONFIGURATION -------------------------------#

# Root directory of the project (edit this to your own path)
PROJECT_ROOT="/path/to/project_root"

# Bowtie2 GRCh38 index prefix (edit this to your own path)
# e.g. if your index files are hg38.1.bt2, hg38.2.bt2, ... then:
# BT2_INDEX_PREFIX="/path/to/bowtie2_indexes/hg38"
BT2_INDEX_PREFIX="/path/to/bowtie2_indexes/GRCh38"

# Number of threads for Bowtie2 and samtools
THREADS=8

# ------------------------- INPUT ARGUMENTS ---------------------------------#

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <COHORT>"
  echo "  COHORT must be one of: CAD, TB"
  exit 1
fi

COHORT="$1"

case "${COHORT}" in
  CAD|cad)
    COHORT="CAD"
    ;;
  TB|tb)
    COHORT="TB"
    ;;
  *)
    echo "Error: COHORT must be CAD or TB"
    exit 1
    ;;
esac

echo ">>> Running host mapping and non-host extraction for cohort: ${COHORT}"

# ------------------------- DIRECTORY LAYOUT --------------------------------#

# Input FASTQ directory (after fastp QC)
# You can adjust subdirectory names if needed.
FASTQ_DIR="${PROJECT_ROOT}/fastq_qc/${COHORT}"

# Output directories
HOST_BAM_DIR="${PROJECT_ROOT}/hostmap/${COHORT}"
UNM_STRICT_DIR="${PROJECT_ROOT}/nonhost/${COHORT}/strict_paired_unmapped"
UNM_RELAX_DIR="${PROJECT_ROOT}/nonhost/${COHORT}/relaxed_any_unmapped"
LOG_DIR="${PROJECT_ROOT}/logs/hostmap_${COHORT}"

mkdir -p "${HOST_BAM_DIR}" "${UNM_STRICT_DIR}" "${UNM_RELAX_DIR}" "${LOG_DIR}"

# ------------------------- MAIN LOOP ---------------------------------------#

for fq1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
  # Skip if no files found
  if [[ ! -e "${fq1}" ]]; then
    echo "No FASTQ files found in ${FASTQ_DIR}"
    exit 1
  fi

  sample=$(basename "${fq1}" _R1.fastq.gz)
  fq2="${FASTQ_DIR}/${sample}_R2.fastq.gz"

  if [[ ! -f "${fq2}" ]]; then
    echo "Warning: R2 file missing for sample ${sample}, skipping."
    continue
  fi

  echo "====================================================================="
  echo ">>> Processing sample: ${sample}"
  echo "====================================================================="

  # -------------------- 1. Bowtie2 mapping to human GRCh38 --------------#
  host_sam="${HOST_BAM_DIR}/${sample}.host.sam"
  host_bam="${HOST_BAM_DIR}/${sample}.host.bam"
  host_bam_sorted="${HOST_BAM_DIR}/${sample}.host.sorted.bam"

  if [[ -f "${host_bam_sorted}" ]]; then
    echo ">>> [${sample}] Sorted host BAM already exists, skipping alignment."
  else
    echo ">>> [${sample}] Bowtie2 mapping to GRCh38"

    bowtie2 \
      -x "${BT2_INDEX_PREFIX}" \
      -1 "${fq1}" \
      -2 "${fq2}" \
      --very-sensitive \
      -p "${THREADS}" \
      -S "${host_sam}" \
      2> "${LOG_DIR}/${sample}.bowtie2.log"

    echo ">>> [${sample}] Convert SAM -> BAM"
    samtools view -@ "${THREADS}" -bS "${host_sam}" > "${host_bam}"

    echo ">>> [${sample}] Sort BAM"
    samtools sort -@ "${THREADS}" -o "${host_bam_sorted}" "${host_bam}"

    echo ">>> [${sample}] Index BAM"
    samtools index "${host_bam_sorted}"

    # Remove intermediate files to save space
    rm -f "${host_sam}" "${host_bam}"
  fi

  # -------------------- 2. STRICT non-host (paired-unmapped) -------------#
  strict_bam="${HOST_BAM_DIR}/${sample}.strict_paired_unmapped.bam"

  if [[ -f "${strict_bam}" ]]; then
    echo ">>> [${sample}] Strict BAM already exists, skipping extraction."
  else
    echo ">>> [${sample}] Extract strictly non-host pairs (both mates unmapped)"

    # -f 12 : read unmapped (0x4) AND mate unmapped (0x8)
    # -F 256: remove secondary alignments
    samtools view -@ "${THREADS}" -b -f 12 -F 256 \
      "${host_bam_sorted}" > "${strict_bam}"
  fi

  echo ">>> [${sample}] Convert strict BAM -> FASTQ"
  samtools fastq \
    -@ "${THREADS}" \
    -1 "${UNM_STRICT_DIR}/${sample}.strict_R1.fastq.gz" \
    -2 "${UNM_STRICT_DIR}/${sample}.strict_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n \
    "${strict_bam}"

  # -------------------- 3. RELAXED non-host (any-unmapped) ---------------#
  relaxed_bam="${HOST_BAM_DIR}/${sample}.relaxed_any_unmapped.bam"

  if [[ -f "${relaxed_bam}" ]]; then
    echo ">>> [${sample}] Relaxed BAM already exists, skipping extraction."
  else
    echo ">>> [${sample}] Extract relaxed non-host reads (any mate unmapped)"

    # "Any-unmapped" definition:
    #   include reads where this read is unmapped (0x4),
    #   regardless of mate status.
    # -f 4  : read unmapped
    # -F 260: remove secondary (0x100) and QC-fail (0x200) reads
    samtools view -@ "${THREADS}" -b -f 4 -F 260 \
      "${host_bam_sorted}" > "${relaxed_bam}"
  fi

  echo ">>> [${sample}] Convert relaxed BAM -> FASTQ"
  samtools fastq \
    -@ "${THREADS}" \
    -1 "${UNM_RELAX_DIR}/${sample}.relaxed_R1.fastq.gz" \
    -2 "${UNM_RELAX_DIR}/${sample}.relaxed_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    -n \
    "${relaxed_bam}"

  echo ">>> [${sample}] DONE."
done

echo ">>> All samples for cohort ${COHORT} finished successfully."
