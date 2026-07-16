#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 05_build_hmms.sh
# hmmbuild from each family's trimmed TRAIN alignment, then score the
# positive test set and hard-negative set against it (used by
# 06_calibrate_cutoffs.py). Assumes hmmer is on PATH.
#
# NOTE: hmmbuild is given -n "${FAMILY}" explicitly so the resulting HMM's
# model name is the family code (e.g. "ectC"), not the alignment file's
# internal MSA name (e.g. "ectC.trimmed") — otherwise downstream family
# matching in compute_metrics.py silently fails (all HMM calls attributed
# to the wrong "family").
#
# Usage: bash 05_build_hmms.sh [refs_dir] [alignments_dir] [hmms_dir]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
ALN_DIR="${2:-alignments}"
HMM_DIR="${3:-hmms}"
THREADS="${THREADS:-4}"

mkdir -p "${HMM_DIR}/scores"

for TRIMMED_ALN in "${ALN_DIR}"/*.trimmed.fasta; do
    [ -s "${TRIMMED_ALN}" ] || continue
    FAMILY=$(basename "${TRIMMED_ALN}" .trimmed.fasta)

    HMM_FILE="${HMM_DIR}/${FAMILY}.hmm"
    POS_SCORES="${HMM_DIR}/scores/${FAMILY}.positive.tblout"
    NEG_SCORES="${HMM_DIR}/scores/${FAMILY}.negative.tblout"

    echo "[${FAMILY}] hmmbuild..."
    hmmbuild --amino -n "${FAMILY}" --cpu "${THREADS}" "${HMM_FILE}" "${TRIMMED_ALN}" > /dev/null

    # Score against the held-out positive TEST set (not train) as a sanity
    # check of true-positive recall on unseen data.
    if [ -s "${REFS_DIR}/${FAMILY}.positive.test.faa" ]; then
        hmmsearch --cpu "${THREADS}" --tblout "${POS_SCORES}" --noali \
            "${HMM_FILE}" "${REFS_DIR}/${FAMILY}.positive.test.faa" > /dev/null
    elif [ -s "${REFS_DIR}/${FAMILY}.positive.faa" ]; then
        hmmsearch --cpu "${THREADS}" --tblout "${POS_SCORES}" --noali \
            "${HMM_FILE}" "${REFS_DIR}/${FAMILY}.positive.faa" > /dev/null
    fi

    if [ -s "${REFS_DIR}/${FAMILY}.negative.faa" ]; then
        hmmsearch --cpu "${THREADS}" --tblout "${NEG_SCORES}" --noali \
            "${HMM_FILE}" "${REFS_DIR}/${FAMILY}.negative.faa" > /dev/null
    else
        : > "${NEG_SCORES}"   # empty; 06_calibrate_cutoffs.py handles this gracefully
    fi

    echo "[${FAMILY}] done -> ${HMM_FILE}"
done

echo ""
echo "hmmbuild + scoring complete. Next: python 06_calibrate_cutoffs.py"
