#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 02_cluster_cdhit.sh
# Cluster each family's positive AND negative FASTA at 90% amino acid
# identity to remove redundancy before the train/test split. Clustering
# negatives too (not just positives) matters for the same anti-leakage
# reason: 03_split_train_test.py splits negative.faa for decoy_from_negatives
# families, and near-duplicate negatives split randomly across train/test
# would let a decoy DIAMOND ref "see" a near-identical sequence to one of its
# own benchmark negatives, inflating apparent decoy performance. Assumes
# cd-hit is on PATH (true inside the osmo_refdb Docker image / conda env).
#
# Usage: bash 02_cluster_cdhit.sh [refs_dir]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
THREADS="${THREADS:-4}"

for LABEL in positive negative; do
    for FASTA in "${REFS_DIR}"/*."${LABEL}".faa; do
        [ -s "${FASTA}" ] || continue
        FAMILY=$(basename "${FASTA}" ".${LABEL}.faa")
        CLUSTERED="${REFS_DIR}/${FAMILY}.${LABEL}.clustered.faa"

        N_SEQ=$(grep -c '^>' "${FASTA}")
        if [ "${N_SEQ}" -lt 5 ]; then
            echo "[${FAMILY}/${LABEL}] SKIP clustering: only ${N_SEQ} sequences, keeping as-is"
            cp "${FASTA}" "${CLUSTERED}"
            continue
        fi

        echo "[${FAMILY}/${LABEL}] clustering ${N_SEQ} sequences at 90% identity..."
        cd-hit -i "${FASTA}" -o "${CLUSTERED}" -c 0.90 -n 5 -T "${THREADS}" -M 16000 -d 0 > /dev/null
        N_CLUSTERED=$(grep -c '^>' "${CLUSTERED}")
        echo "[${FAMILY}/${LABEL}] ${N_SEQ} -> ${N_CLUSTERED} sequences after clustering"

        # Overwrite <label>.faa with the deduplicated set so downstream steps
        # (03_split_train_test.py onward) only ever see non-redundant sequences.
        cp "${CLUSTERED}" "${FASTA}"
    done
done

echo ""
echo "Clustering complete."
