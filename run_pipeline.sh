#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — run_pipeline.sh
#
# Single entrypoint for building and benchmarking osmo_refdb reference
# databases (DIAMOND + HMM) for osmoadaptation gene families.
#
# Usage:
#   bash run_pipeline.sh build [release_name]      # fetch -> cluster -> split
#                                                   # -> align -> hmmbuild ->
#                                                   # calibrate -> press ->
#                                                   # diamond makedb
#   bash run_pipeline.sh benchmark [release_name]  # simulate reads -> run
#                                                   # DIAMOND+HMM -> metrics
#   bash run_pipeline.sh all [release_name]        # build then benchmark
#
# All outputs for a given release_name (default: "dev") are written under
# releases/<release_name>/ so multiple database versions can coexist and be
# reproduced independently (e.g. "v1", "v2-added-ectD").
#
# Adding a new gene family: edit families.yaml, then rerun `build`.
#
# Run this from inside the osmo_refdb Docker container or a conda env with
# environment.yml activated.
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

CMD="${1:-}"
RELEASE_NAME="${2:-dev}"
RELEASE_DIR="releases/${RELEASE_NAME}"
THREADS="${THREADS:-4}"

REFS_DIR="${RELEASE_DIR}/refs"
ALN_DIR="${RELEASE_DIR}/alignments"
HMM_DIR="${RELEASE_DIR}/hmms"
RESULTS_DIR="${RELEASE_DIR}/results"

usage() {
    echo "Usage: bash run_pipeline.sh {build|benchmark|all} [release_name]"
    exit 1
}

do_build() {
    echo "### osmo_refdb build — release '${RELEASE_NAME}' ###"
    mkdir -p "${REFS_DIR}" "${ALN_DIR}" "${HMM_DIR}"

    echo ""
    echo "--- 1. Fetch reference sequences from UniProt ---"
    python pipeline/01_fetch_refs.py --families families.yaml --out "${REFS_DIR}"

    echo ""
    echo "--- 2. CD-HIT cluster positives (remove redundancy) ---"
    bash pipeline/02_cluster_cdhit.sh "${REFS_DIR}"

    echo ""
    echo "--- 3. Train/test split (held-out benchmark set, no leakage) ---"
    python pipeline/03_split_train_test.py --refs "${REFS_DIR}" --families families.yaml

    echo ""
    echo "--- 4. Align + trim (TRAIN positives only) ---"
    THREADS="${THREADS}" bash pipeline/04_align_trim.sh "${REFS_DIR}" "${ALN_DIR}"

    echo ""
    echo "--- 5. hmmbuild + score positive-test/negative sets ---"
    THREADS="${THREADS}" bash pipeline/05_build_hmms.sh "${REFS_DIR}" "${ALN_DIR}" "${HMM_DIR}"

    echo ""
    echo "--- 6. Calibrate GA cutoffs ---"
    python pipeline/06_calibrate_cutoffs.py --hmms "${HMM_DIR}" --families families.yaml

    echo ""
    echo "--- 7. Press combined HMM database ---"
    bash pipeline/07_press_hmms.sh "${HMM_DIR}" osmo_refdb

    echo ""
    echo "--- 8. Build DIAMOND database (TRAIN positives only) ---"
    THREADS="${THREADS}" bash pipeline/08_build_diamond_db.sh "${REFS_DIR}" "${RELEASE_DIR}" osmo_refdb

    echo ""
    echo "=== Build complete: ${RELEASE_DIR} ==="
    echo "  DIAMOND db : ${RELEASE_DIR}/osmo_refdb.dmnd"
    echo "  HMM db     : ${HMM_DIR}/osmo_refdb.hmm"
    echo "  Cutoffs    : ${HMM_DIR}/cutoff_manifest.tsv"
    echo "Next: bash run_pipeline.sh benchmark ${RELEASE_NAME}"
}

do_benchmark() {
    echo "### osmo_refdb benchmark — release '${RELEASE_NAME}' ###"
    if [ ! -s "${RELEASE_DIR}/osmo_refdb.dmnd" ]; then
        echo "ERROR: ${RELEASE_DIR}/osmo_refdb.dmnd not found. Run 'build' first." >&2
        exit 1
    fi
    if [ ! -s "${HMM_DIR}/osmo_refdb.hmm" ]; then
        echo "ERROR: ${HMM_DIR}/osmo_refdb.hmm not found. Run 'build' first." >&2
        exit 1
    fi

    mkdir -p "${RESULTS_DIR}"

    echo ""
    echo "--- 9. Simulate reads from held-out test positives + hard negatives ---"
    python pipeline/09_simulate_reads.py \
        --refs "${REFS_DIR}" --families families.yaml \
        --reads-per-sequence 20 \
        --out "${RESULTS_DIR}/reads"

    echo ""
    echo "--- 10. Run DIAMOND (train-only db) vs HMM (train-only, calibrated) ---"
    THREADS="${THREADS}" bash pipeline/10_run_benchmark.sh \
        "${RESULTS_DIR}/reads" "${RESULTS_DIR}" \
        "${RELEASE_DIR}/osmo_refdb.dmnd" "${HMM_DIR}/osmo_refdb.hmm"

    echo ""
    echo "--- 11. Compute precision/recall/F1 + ROC/PR curves ---"
    python pipeline/11_compute_metrics.py \
        --results "${RESULTS_DIR}" --out "${RESULTS_DIR}/metrics" --families families.yaml

    echo ""
    echo "=== Benchmark complete. See ${RESULTS_DIR}/metrics for summary.tsv, "
    echo "    best_threshold_summary.tsv, and ROC/PR plots. ==="
}

case "${CMD}" in
    build)
        do_build
        ;;
    benchmark)
        do_benchmark
        ;;
    all)
        do_build
        do_benchmark
        ;;
    *)
        usage
        ;;
esac
