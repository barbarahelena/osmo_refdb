#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 10_run_benchmark.sh
# Run DIAMOND (via osmotool profile, train-only db) and HMM (translate +
# hmmscan, train-only calibrated HMMs) on the simulated reads, for later
# scoring by 11_compute_metrics.py.
#
# Usage: bash 10_run_benchmark.sh [reads_dir] [results_dir] [diamond_db] [hmm_db]
# ---------------------------------------------------------------------------

set -euo pipefail

READS_DIR="${1:-results/reads}"
OUT_DIR="${2:-results}"
DIAMOND_DB="${3:-releases/dev/osmo_refdb.dmnd}"
HMM_DB="${4:-hmms/osmo_refdb.hmm}"
THREADS="${THREADS:-4}"
# The benchmark's positive test set is deliberately broad/phylogenetically
# diverse (all UniProt hits for the gene, not just close orthologs), so
# many true positives align at only 50-70% identity to their nearest
# training-set reference. osmotool profile's published default (--id 80)
# is tuned for high-confidence calls in typical use, but would silently
# drop most divergent true positives here. Relax to more permissive
# defaults for benchmarking; override via env vars if needed.
DIAMOND_MIN_IDENTITY="${DIAMOND_MIN_IDENTITY:-0.40}"
DIAMOND_MIN_QUERY_COVER="${DIAMOND_MIN_QUERY_COVER:-0.50}"
# Off by default (0.0), same as osmotool's own default -- set this to test
# whether filtering on subject/target coverage reduces the multi-domain
# fusion-protein false positives found while benchmarking proX (see
# --min_subject_cover's help text in osmotool for why this isn't on by
# default: a short read is always much shorter than most full-length
# reference proteins, so a strict threshold here also costs real recall).
DIAMOND_MIN_SUBJECT_COVER="${DIAMOND_MIN_SUBJECT_COVER:-0.0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${OUT_DIR}/diamond" "${OUT_DIR}/hmm"

echo "=== DIAMOND (osmotool profile, --min_identity ${DIAMOND_MIN_IDENTITY} --min_query_cover ${DIAMOND_MIN_QUERY_COVER} --min_subject_cover ${DIAMOND_MIN_SUBJECT_COVER}) ==="
for R1 in "${READS_DIR}"/*_R1.fastq*; do
    [ -e "${R1}" ] || continue
    SAMPLE=$(basename "${R1}" | sed -E 's/_R1\.fastq(\.gz)?$//')
    R2="${READS_DIR}/${SAMPLE}_R2.fastq"
    [ -f "${R2}.gz" ] && R2="${R2}.gz"
    [ -f "${R2}" ] || { echo "  SKIP ${SAMPLE}: no R2 found"; continue; }

    echo "  sample: ${SAMPLE}"
    osmotool profile \
        "${DIAMOND_DB}" \
        -1 "${R1}" -2 "${R2}" \
        --out_prefix "${OUT_DIR}/diamond/${SAMPLE}" \
        --threads "${THREADS}" \
        --min_identity "${DIAMOND_MIN_IDENTITY}" \
        --min_query_cover "${DIAMOND_MIN_QUERY_COVER}" \
        --min_subject_cover "${DIAMOND_MIN_SUBJECT_COVER}" \
        --keep_aln
done

echo ""
echo "=== HMM (6-frame translate + hmmscan) ==="
for R1 in "${READS_DIR}"/*_R1.fastq*; do
    [ -e "${R1}" ] || continue
    SAMPLE=$(basename "${R1}" | sed -E 's/_R1\.fastq(\.gz)?$//')
    R2="${READS_DIR}/${SAMPLE}_R2.fastq"
    [ -f "${R2}.gz" ] && R2="${R2}.gz"
    [ -f "${R2}" ] || { echo "  SKIP ${SAMPLE}: no R2 found"; continue; }

    echo "  sample: ${SAMPLE}"
    bash "${SCRIPT_DIR}/translate_and_hmmscan.sh" \
        "${R1}" "${R2}" "${HMM_DB}" "${OUT_DIR}/hmm/${SAMPLE}" "${THREADS}"
done

echo ""
echo "Done. Now run: python 11_compute_metrics.py --results ${OUT_DIR} --out ${OUT_DIR}/metrics"
