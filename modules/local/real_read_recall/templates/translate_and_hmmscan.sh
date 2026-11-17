#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — helper: translate reads (orfm 6-frame) + hmmscan
# Usage: translate_and_hmmscan.sh R1.fastq R2.fastq osmo_refdb.hmm out_prefix threads
# ---------------------------------------------------------------------------

set -euo pipefail

R1="$1"
R2="$2"
HMM_DB="$3"
OUT_PREFIX="$4"
THREADS="${5:-4}"

MERGED_FASTQ="${OUT_PREFIX}.merged.fastq"
ORFS_FASTA="${OUT_PREFIX}.orfs.faa"
TBLOUT="${OUT_PREFIX}.hmmscan.tblout"

mkdir -p "$(dirname "${OUT_PREFIX}")"

zcat -f "${R1}" "${R2}" > "${MERGED_FASTQ}"

orfm "${MERGED_FASTQ}" > "${ORFS_FASTA}"

# NOTE: no --cut_ga here. GA cutoffs were calibrated on full-length proteins
# (see 06_calibrate_cutoffs.py); short simulated-read ORF fragments can never
# reach that bit-score threshold, so using --cut_ga would silently produce
# zero hits for every read. Instead emit raw hmmscan bit scores and let
# 11_compute_metrics.py sweep thresholds (same approach used for DIAMOND).
hmmscan --cpu "${THREADS}" --tblout "${TBLOUT}" --noali \
    "${HMM_DB}" "${ORFS_FASTA}" > /dev/null

rm -f "${MERGED_FASTQ}"
echo "  HMM hits: ${TBLOUT}"
