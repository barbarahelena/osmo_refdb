#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 04_align_trim.sh
# MAFFT alignment + trimAl trimming of each family's TRAIN positives.
# Assumes mafft and trimal are on PATH (true inside the osmo_refdb Docker
# image / conda env, or an HPC module-loaded shell — no SLURM/Singularity
# wrapping needed since this step is fast even run as a simple loop).
#
# trimAl's -gt (gap threshold) defaults to 0.8 for every family, but can be
# overridden per family via families.yaml's `trim_gt` field (see ectA/ectB:
# an experiment testing whether their broad, promiscuous Pfam folds need
# more of the alignment's flanking columns kept, not trimmed away, to
# preserve their true discriminating signal).
#
# Usage: bash 04_align_trim.sh [refs_dir] [alignments_dir] [families_yaml]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
ALN_DIR="${2:-alignments}"
FAMILIES_FILE="${3:-families.yaml}"
THREADS="${THREADS:-4}"

mkdir -p "${ALN_DIR}"

trim_gt_for_family() {
    python3 -c "
import sys, yaml
with open(sys.argv[1]) as fh:
    data = yaml.safe_load(fh)
for fam in data['families']:
    if fam['name'] == sys.argv[2]:
        print(fam.get('trim_gt', 0.8))
        break
else:
    print(0.8)
" "${FAMILIES_FILE}" "$1"
}

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

    GT=$(trim_gt_for_family "${FAMILY}")
    mafft --auto --thread "${THREADS}" "${INPUT}" > "${RAW_ALN}" 2>/dev/null
    trimal -in "${RAW_ALN}" -out "${TRIMMED_ALN}" -gt "${GT}"

    echo "[${FAMILY}] done (trim_gt=${GT}) -> ${TRIMMED_ALN}"
done

echo ""
echo "Alignment step complete."
