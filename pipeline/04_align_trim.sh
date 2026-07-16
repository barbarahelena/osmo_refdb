#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 04_align_trim.sh
# MAFFT alignment + trimAl trimming of each family's TRAIN positives.
# Assumes mafft and trimal are on PATH (true inside the osmo_refdb Docker
# image / conda env, or an HPC module-loaded shell — no SLURM/Singularity
# wrapping needed since this step is fast even run as a simple loop).
#
# Usage: bash 04_align_trim.sh [refs_dir] [alignments_dir]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
ALN_DIR="${2:-alignments}"
THREADS="${THREADS:-4}"

mkdir -p "${ALN_DIR}"

for INPUT in "${REFS_DIR}"/*.positive.train.faa; do
    [ -s "${INPUT}" ] || continue
    FAMILY=$(basename "${INPUT}" .positive.train.faa)

    RAW_ALN="${ALN_DIR}/${FAMILY}.aln.fasta"
    TRIMMED_ALN="${ALN_DIR}/${FAMILY}.trimmed.fasta"

    N_SEQ=$(grep -c '^>' "${INPUT}")
    echo "[${FAMILY}] aligning ${N_SEQ} sequences..."

    if [ "${N_SEQ}" -lt 2 ]; then
        echo "[${FAMILY}] SKIP: need >=2 sequences to align, found ${N_SEQ}"
        continue
    fi

    mafft --auto --thread "${THREADS}" "${INPUT}" > "${RAW_ALN}" 2>/dev/null
    trimal -in "${RAW_ALN}" -out "${TRIMMED_ALN}" -gt 0.8

    echo "[${FAMILY}] done -> ${TRIMMED_ALN}"
done

echo ""
echo "Alignment step complete."
