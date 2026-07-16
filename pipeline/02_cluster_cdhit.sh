#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 02_cluster_cdhit.sh
# Cluster each family's positive FASTA at 90% amino acid identity to remove
# redundancy before alignment/HMM building. Assumes cd-hit is on PATH
# (true inside the osmo_refdb Docker image / conda env).
#
# Usage: bash 02_cluster_cdhit.sh [refs_dir]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
THREADS="${THREADS:-4}"

for POS_FASTA in "${REFS_DIR}"/*.positive.faa; do
    [ -s "${POS_FASTA}" ] || continue
    FAMILY=$(basename "${POS_FASTA}" .positive.faa)
    CLUSTERED="${REFS_DIR}/${FAMILY}.positive.clustered.faa"

    N_SEQ=$(grep -c '^>' "${POS_FASTA}")
    if [ "${N_SEQ}" -lt 5 ]; then
        echo "[${FAMILY}] SKIP clustering: only ${N_SEQ} sequences, keeping as-is"
        cp "${POS_FASTA}" "${CLUSTERED}"
        continue
    fi

    echo "[${FAMILY}] clustering ${N_SEQ} sequences at 90% identity..."
    cd-hit -i "${POS_FASTA}" -o "${CLUSTERED}" -c 0.90 -n 5 -T "${THREADS}" -M 16000 -d 0 > /dev/null
    N_CLUSTERED=$(grep -c '^>' "${CLUSTERED}")
    echo "[${FAMILY}] ${N_SEQ} -> ${N_CLUSTERED} sequences after clustering"

    # Overwrite positive.faa with the deduplicated set so downstream steps
    # (03_split_train_test.py onward) only ever see non-redundant sequences.
    cp "${CLUSTERED}" "${POS_FASTA}"
done

echo ""
echo "Clustering complete."
