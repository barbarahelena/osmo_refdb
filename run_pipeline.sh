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
# SPLIT_MODE=random|taxonomy (default random) controls how train/test is
# split in step 3 -- "taxonomy" holds out whole genera instead of individual
# sequences, testing detection of genuinely divergent homologs rather than
# just held-out sequences from the same overall pool.
#
# FAMILIES_FILE=families.yaml (default) lets you point the whole pipeline
# at a subset file instead -- e.g. to test-drive one or two newly-added
# families (fetch -> QC -> calibrate -> benchmark, ~minutes) before
# committing to a full rebuild of every family (~90 minutes). Use a
# DIFFERENT release_name for the subset run (e.g. "v3-test") so it writes
# to its own releases/ directory and can't leave the real release in a
# partial state -- see pipeline/make_family_subset.py to generate the
# subset file, and the README for the full workflow.
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
FAMILIES_FILE="${FAMILIES_FILE:-families.yaml}"

REFS_DIR="${RELEASE_DIR}/refs"
ALN_DIR="${RELEASE_DIR}/alignments"
HMM_DIR="${RELEASE_DIR}/hmms"
RESULTS_DIR="${RELEASE_DIR}/results"
SPLIT_MODE="${SPLIT_MODE:-random}"

usage() {
    echo "Usage: bash run_pipeline.sh {build|benchmark|all} [release_name]"
    exit 1
}

do_build() {
    echo "### osmo_refdb build — release '${RELEASE_NAME}' ###"
    mkdir -p "${REFS_DIR}" "${ALN_DIR}" "${HMM_DIR}"

    echo ""
    echo "--- 1. Fetch reference sequences from UniProt ---"
    python pipeline/01_fetch_refs.py --families "${FAMILIES_FILE}" --out "${REFS_DIR}"

    echo ""
    echo "--- 1b. Check hard-negative purity (flag negatives that are secretly true positives) ---"
    python pipeline/01b_check_negative_purity.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}" --threads "${THREADS}"

    echo ""
    echo "--- 1c. Check length outliers (flag likely fusion proteins/fragments) ---"
    python pipeline/01c_check_length_outliers.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}"

    echo ""
    echo "--- 1d. Merge extra positives from local studies (if any) ---"
    python pipeline/01d_add_extra_positives.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}"

    echo ""
    echo "--- 2. CD-HIT cluster positives + negatives (remove redundancy) ---"
    bash pipeline/02_cluster_cdhit.sh "${REFS_DIR}"

    echo ""
    echo "--- 3. Train/test split (held-out benchmark set, no leakage; mode=${SPLIT_MODE}) ---"
    python pipeline/03_split_train_test.py --refs "${REFS_DIR}" --families "${FAMILIES_FILE}" \
        --split-mode "${SPLIT_MODE}"

    echo ""
    echo "--- 4. Align + trim (TRAIN positives only) ---"
    THREADS="${THREADS}" bash pipeline/04_align_trim.sh "${REFS_DIR}" "${ALN_DIR}" "${FAMILIES_FILE}"

    echo ""
    echo "--- 5. hmmbuild (or fetch curated Pfam model) + score positive-test/negative sets ---"
    THREADS="${THREADS}" bash pipeline/05_build_hmms.sh "${REFS_DIR}" "${ALN_DIR}" "${HMM_DIR}" "${FAMILIES_FILE}"

    echo ""
    echo "--- 6. Calibrate GA cutoffs ---"
    python pipeline/06_calibrate_cutoffs.py --hmms "${HMM_DIR}" --families "${FAMILIES_FILE}"

    echo ""
    echo "--- 6b. Consolidated QC scorecard ---"
    python pipeline/06b_qc_scorecard.py --refs "${REFS_DIR}" --hmms "${HMM_DIR}" \
        --families "${FAMILIES_FILE}" --out "${RELEASE_DIR}/qc_scorecard.tsv"

    echo ""
    echo "--- 7. Press combined HMM database ---"
    bash pipeline/07_press_hmms.sh "${HMM_DIR}" osmo_refdb

    echo ""
    echo "--- 8a. Build decoy references (confusable paralogs, e.g. betL vs betT/caiT) ---"
    python pipeline/08a_build_decoy_refs.py --refs "${REFS_DIR}" --families "${FAMILIES_FILE}"

    echo ""
    echo "--- 8d. Build fused-ORF references (e.g. mrpA/mrpB single-ORF lineages) ---"
    python pipeline/08d_build_fusion_refs.py --refs "${REFS_DIR}" --families "${FAMILIES_FILE}"

    echo ""
    echo "--- 8. Build DIAMOND database (TRAIN positives + any decoy/fused refs) ---"
    THREADS="${THREADS}" bash pipeline/08_build_diamond_db.sh "${REFS_DIR}" "${RELEASE_DIR}" osmo_refdb

    echo ""
    echo "--- 8b. Calibrate per-family DIAMOND bit-score cutoffs ---"
    python pipeline/08b_calibrate_diamond_cutoffs.py \
        --refs "${REFS_DIR}" --release "${RELEASE_DIR}" --release-name osmo_refdb \
        --families "${FAMILIES_FILE}" --threads "${THREADS}"

    echo ""
    echo "--- 8c. Write scope exclusion lists (annotate_only + decoy families) ---"
    python pipeline/08c_write_scope_manifest.py --families "${FAMILIES_FILE}" \
        --profile-out "${RELEASE_DIR}/osmo_refdb.profile_excluded_families.txt" \
        --annotate-out "${RELEASE_DIR}/osmo_refdb.annotate_excluded_families.txt"

    echo ""
    echo "=== Build complete: ${RELEASE_DIR} ==="
    echo "  DIAMOND db : ${RELEASE_DIR}/osmo_refdb.dmnd"
    echo "  HMM db     : ${HMM_DIR}/osmo_refdb.hmm"
    echo "  Cutoffs    : ${HMM_DIR}/cutoff_manifest.tsv"
    echo "  DIAMOND cutoffs: ${RELEASE_DIR}/osmo_refdb.diamond_cutoffs.tsv"
    echo "  QC scorecard: ${RELEASE_DIR}/qc_scorecard.tsv"
    echo "  Profile-excluded families: ${RELEASE_DIR}/osmo_refdb.profile_excluded_families.txt"
    echo "  Annotate-excluded families: ${RELEASE_DIR}/osmo_refdb.annotate_excluded_families.txt"
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
    # Optional: READ_LENGTHS="100,150,250,300" for read-length-stratified
    # metrics; BACKGROUND_R1/BACKGROUND_R2 to score a real metagenome
    # sample's false-positive rate (see background_fpr.tsv).
    EXTRA_SIMULATE_ARGS=()
    [ -n "${READ_LENGTHS:-}" ] && EXTRA_SIMULATE_ARGS+=(--read-lengths "${READ_LENGTHS}")
    [ -n "${BACKGROUND_R1:-}" ] && [ -n "${BACKGROUND_R2:-}" ] && \
        EXTRA_SIMULATE_ARGS+=(--background "${BACKGROUND_R1}" "${BACKGROUND_R2}")
    python pipeline/09_simulate_reads.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}" \
        --reads-per-sequence 20 \
        --out "${RESULTS_DIR}/reads" \
        "${EXTRA_SIMULATE_ARGS[@]+"${EXTRA_SIMULATE_ARGS[@]}"}"

    echo ""
    echo "--- 9b. Compute per-read truth via actual CDS overlap ---"
    python pipeline/09b_compute_read_truth.py --reads "${RESULTS_DIR}/reads"

    echo ""
    echo "--- 10. Run DIAMOND (train-only db) vs HMM (train-only, calibrated) ---"
    THREADS="${THREADS}" bash pipeline/10_run_benchmark.sh \
        "${RESULTS_DIR}/reads" "${RESULTS_DIR}" \
        "${RELEASE_DIR}/osmo_refdb.dmnd" "${HMM_DIR}/osmo_refdb.hmm"

    echo ""
    echo "--- 11. Compute precision/recall/F1 + ROC/PR curves ---"
    python pipeline/11_compute_metrics.py \
        --results "${RESULTS_DIR}" --out "${RESULTS_DIR}/metrics" --families "${FAMILIES_FILE}" \
        --cascade-out "${RELEASE_DIR}/osmo_refdb.profile_cascade.tsv"

    echo ""
    echo "=== Benchmark complete. See ${RESULTS_DIR}/metrics for summary.tsv, "
    echo "    best_threshold_summary.tsv, and ROC/PR plots. ==="
    echo "    Profile-mode DIAMOND+HMM cascade config: ${RELEASE_DIR}/osmo_refdb.profile_cascade.tsv"
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
