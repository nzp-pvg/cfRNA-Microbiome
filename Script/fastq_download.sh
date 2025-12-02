#!/usr/bin/env bash
# -------------------------------------------------------------------
# Script:        01_download_raw_fastq.sh
# Purpose:       Download raw FASTQ files from NCBI SRA and perform
#                basic QC using fastp.
#
# Author:        Zhaoxia Wang
# Version:       v0.9.0
# Date:          2024-09-15
#
# Requirements:
#   - SRA Toolkit  (prefetch, fasterq-dump)
#   - fastp
#   - gzip
#
# Usage:
#   1) Prepare an accession list file, e.g. SRR_Acc_List.txt
#      Each line contains one SRR ID:
#         SRR1234567
#         SRR1234568
#         ...
#   2) Edit the variables below (PROJECT_ROOT, ACC_LIST).
#   3) Run:
#         bash 01_download_raw_fastq.sh
# -------------------------------------------------------------------

set -euo pipefail

# --------- User settings ------------------------------------------------------

# Root directory for this project
PROJECT_ROOT="/path/to/project"        # e.g. /Users/xxx/Science/CVD_MS_5/data

# File containing SRR accessions (one per line)
ACC_LIST="${PROJECT_ROOT}/SRR_Acc_List.txt"

# Directory where SRA Toolkit will store *.sra files
SRA_DIR="${PROJECT_ROOT}/sra"

# Directory for raw FASTQ (after fasterq-dump + gzip)
FASTQ_DIR="${PROJECT_ROOT}/fastq"

# Directory for QC'ed FASTQ (fastp output)
FASTQ_QC_DIR="${PROJECT_ROOT}/fastq_qc"

# Directories for fastp reports and logs
FASTP_HTML_DIR="${PROJECT_ROOT}/fastp_html"
FASTP_JSON_DIR="${PROJECT_ROOT}/fastp_json"
LOG_DIR="${PROJECT_ROOT}/logs"

# Number of CPU threads for fasterq-dump and fastp
N_THREADS=8

# -------------------------------------------------------------------
# Helper: check that a command exists
# -------------------------------------------------------------------
check_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] Command '$1' not found in PATH. Please install it first." >&2
    exit 1
  fi
}

echo ">>> [STEP 0] Checking required tools ..."

check_cmd prefetch
check_cmd fasterq-dump
check_cmd fastp
check_cmd gzip

# -------------------------------------------------------------------
# Create output directories
# -------------------------------------------------------------------
echo ">>> [STEP 1] Creating directories ..."

mkdir -p "${SRA_DIR}" \
         "${FASTQ_DIR}" \
         "${FASTQ_QC_DIR}" \
         "${FASTP_HTML_DIR}" \
         "${FASTP_JSON_DIR}" \
         "${LOG_DIR}"

# -------------------------------------------------------------------
# Download .sra files using prefetch
# -------------------------------------------------------------------
echo ">>> [STEP 2] Downloading SRA files with prefetch ..."
echo "     Accession list: ${ACC_LIST}"
echo "     SRA directory : ${SRA_DIR}"

# prefetch can read the accession list directly via --option-file
prefetch --option-file "${ACC_LIST}" --output-directory "${SRA_DIR}"

echo ">>> [STEP 2] Download completed."

# -------------------------------------------------------------------
# Convert SRA to gzipped paired-end FASTQ with fasterq-dump
# -------------------------------------------------------------------
echo ">>> [STEP 3] Converting SRA to FASTQ (fasterq-dump + gzip) ..."

while read -r SRR; do
  [[ -z "${SRR}" ]] && continue  # skip empty lines

  echo ">>>   Processing ${SRR}"

  fq1="${FASTQ_DIR}/${SRR}_R1.fastq.gz"
  fq2="${FASTQ_DIR}/${SRR}_R2.fastq.gz"

  # Skip if both gzipped FASTQ files already exist
  if [[ -f "${fq1}" && -f "${fq2}" ]]; then
    echo "      FASTQ files already exist, skipping: ${SRR}"
    continue
  fi

  # Temporary uncompressed FASTQ files produced by fasterq-dump
  tmp_fq1="${FASTQ_DIR}/${SRR}_R1.fastq"
  tmp_fq2="${FASTQ_DIR}/${SRR}_R2.fastq"

  # Convert .sra to FASTQ
  fasterq-dump "${SRA_DIR}/${SRR}.sra" \
    --split-files \
    --threads "${N_THREADS}" \
    --outdir "${FASTQ_DIR}"

  # Compress FASTQ files
  gzip -f "${tmp_fq1}"
  gzip -f "${tmp_fq2}"

  echo "      Done: ${SRR}"
done < "${ACC_LIST}"

echo ">>> [STEP 3] FASTQ generation finished."

# -------------------------------------------------------------------
# Run fastp QC on raw FASTQ files
# -------------------------------------------------------------------
echo ">>> [STEP 4] Running fastp quality control ..."

while read -r SRR; do
  [[ -z "${SRR}" ]] && continue  # skip empty lines

  echo ">>>   fastp on ${SRR}"

  raw_fq1="${FASTQ_DIR}/${SRR}_R1.fastq.gz"
  raw_fq2="${FASTQ_DIR}/${SRR}_R2.fastq.gz"

  qc_fq1="${FASTQ_QC_DIR}/${SRR}_R1.fastq.gz"
  qc_fq2="${FASTQ_QC_DIR}/${SRR}_R2.fastq.gz"

  # Skip if QC'ed FASTQ files already exist
  if [[ -f "${qc_fq1}" && -f "${qc_fq2}" ]]; then
    echo "      QC FASTQ already exists, skipping: ${SRR}"
    continue
  fi

  fastp \
    -i "${raw_fq1}" \
    -I "${raw_fq2}" \
    -o "${qc_fq1}" \
    -O "${qc_fq2}" \
    --detect_adapter_for_pe \
    --thread "${N_THREADS}" \
    --html "${FASTP_HTML_DIR}/${SRR}.fastp.html" \
    --json "${FASTP_JSON_DIR}/${SRR}.fastp.json" \
    > "${LOG_DIR}/${SRR}.fastp.log" 2>&1

  echo "      fastp done: ${SRR}"
done < "${ACC_LIST}"

echo ">>> [STEP 4] fastp QC finished."
echo ">>> ALL DONE."