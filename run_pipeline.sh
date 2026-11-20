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
#   bash run_pipeline.sh markers [release_name]    # build + benchmark the
#                                                   # universal single-copy
#                                                   # marker genes
#                                                   # (markers.yaml) as their
#                                                   # own release,
#                                                   # <release_name>-markers,
#                                                   # then package them into
#                                                   # <release_name>/markers/
#                                                   # for osmotool profile's
#                                                   # copies_per_genome column
#   bash run_pipeline.sh catalog [release_name]    # build a gene catalog
#                                                   # (CATALOG_NAME,
#                                                   # CATALOG_GENES) into
#                                                   # <release>/catalogs/ for
#                                                   # osmotool profile --catalog
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
    echo "Usage: bash run_pipeline.sh {build|benchmark|all|markers|nacl_model|catalog} [release_name]"
    exit 1
}

do_markers() {
    # Markers are a separate build on purpose (own DIAMOND db, so they never
    # compete with osmo families for best hits) -- see markers.yaml's header.
    local marker_release="${RELEASE_NAME}-markers"
    if [ ! -s "${RELEASE_DIR}/osmo_refdb.dmnd" ]; then
        echo "ERROR: ${RELEASE_DIR}/osmo_refdb.dmnd not found. Build the main release first." >&2
        exit 1
    fi
    echo "### osmo_refdb markers — building '${marker_release}', packaging into '${RELEASE_NAME}' ###"
    FAMILIES_FILE=markers.yaml SPLIT_MODE="${SPLIT_MODE}" THREADS="${THREADS}" \
        bash "${SCRIPT_DIR}/run_pipeline.sh" all "${marker_release}"
    echo ""
    echo "--- 12. Package markers into ${RELEASE_DIR}/markers/ ---"
    python pipeline/12_package_markers.py \
        --marker-release "releases/${marker_release}" --release "${RELEASE_DIR}" \
        --markers markers.yaml
}

do_nacl_model() {
    # NaCl-max predictor, trained on the BacDive phenotype_validation set
    # annotated with this release (pipeline/phenotype_validation/README.md).
    # Needs xgboost, scikit-learn, pandas, scipy, matplotlib; set PYTHON to an
    # interpreter that has them (see pipeline/nacl_model/README.md).
    local py="${PYTHON:-python}"
    local pv="${RELEASE_DIR}/phenotype_validation"
    local build="${RELEASE_DIR}/nacl_model_build"
    if [ ! -s "${pv}/feature_matrix.tsv" ]; then
        echo "ERROR: ${pv}/feature_matrix.tsv not found. Run the phenotype validation first." >&2
        exit 1
    fi
    echo "### osmo_refdb NaCl-max model — release '${RELEASE_NAME}' ###"
    echo "--- 13a. Train (nested genus/family-blocked CV) ---"
    "${py}" pipeline/nacl_model/train.py --pv-dir "${pv}" --release "${RELEASE_DIR}" \
        --out "${build}" --threads "${THREADS}"
    echo "--- 13b. Explain ---"
    "${py}" pipeline/nacl_model/explain.py --build "${build}" --pv-dir "${pv}"
    echo "--- 13c. Package into ${RELEASE_DIR}/nacl_model/ ---"
    "${py}" pipeline/13_package_nacl_model.py --build "${build}" --release "${RELEASE_DIR}"
}

do_catalog() {
    # Gene catalog for profile --catalog, built into the release with this
    # release's HMMs and shipped with it as one unit (pipeline/14_package_catalog.py).
    # CATALOG_NAME (e.g. gmgc_gut) and CATALOG_GENES (nucleotide CDS FASTA)
    # are required; CATALOG_PROTEINS, CATALOG_SOURCE, CATALOG_SOURCE_URL and
    # CATALOG_INDEX (hits|all|none, default hits: small enough for the
    # release; the full index is built separately with `osmotool
    # index-catalog`) are optional. Needs osmotool with build-catalog and
    # bwa-mem3 on PATH.
    if [ -z "${CATALOG_NAME:-}" ] || [ -z "${CATALOG_GENES:-}" ]; then
        echo "ERROR: set CATALOG_NAME and CATALOG_GENES (see do_catalog in run_pipeline.sh)." >&2
        exit 1
    fi
    if [ ! -s "${HMM_DIR}/osmo_refdb.hmm" ]; then
        echo "ERROR: ${HMM_DIR}/osmo_refdb.hmm not found. Build the release first." >&2
        exit 1
    fi
    echo "### osmo_refdb catalog '${CATALOG_NAME}' — release '${RELEASE_NAME}' ###"
    echo "--- 14. Build + package ${RELEASE_DIR}/catalogs/${CATALOG_NAME}/ ---"
    EXTRA_CATALOG_ARGS=()
    [ -n "${CATALOG_PROTEINS:-}" ] && EXTRA_CATALOG_ARGS+=(--proteins "${CATALOG_PROTEINS}")
    [ "${CATALOG_FORCE:-0}" = "1" ] && EXTRA_CATALOG_ARGS+=(--force)
    python pipeline/14_package_catalog.py \
        --release "${RELEASE_DIR}" --name "${CATALOG_NAME}" --genes "${CATALOG_GENES}" \
        --source "${CATALOG_SOURCE:-}" --source-url "${CATALOG_SOURCE_URL:-}" \
        --index "${CATALOG_INDEX:-hits}" --threads "${THREADS}" \
        "${EXTRA_CATALOG_ARGS[@]+"${EXTRA_CATALOG_ARGS[@]}"}"
}

do_build() {
    echo "### osmo_refdb build — release '${RELEASE_NAME}' ###"
    mkdir -p "${REFS_DIR}" "${ALN_DIR}" "${HMM_DIR}"

    if [ -n "${REUSE_REFS_FROM:-}" ]; then
        # Rebuild from another release's fetched + QC'd + clustered
        # sequences (steps 1-2), so that only the split, calibration and
        # later steps differ and no UniProt/RefSeq drift enters. Copies the
        # whole refs/ dir; step 3 onward overwrites every derived file.
        echo ""
        echo "--- 1-2. Reusing fetched and clustered refs from ${REUSE_REFS_FROM} ---"
        [ -d "${REUSE_REFS_FROM}/refs" ] || { echo "ERROR: ${REUSE_REFS_FROM}/refs not found" >&2; exit 1; }
        cp -R "${REUSE_REFS_FROM}/refs/." "${REFS_DIR}/"
        rm -f "${REFS_DIR}"/*.positive.train.faa "${REFS_DIR}"/*.positive.calib.faa \
              "${REFS_DIR}"/*.positive.test.faa "${REFS_DIR}"/*.negative.train.faa \
              "${REFS_DIR}"/*.negative.test.faa "${REFS_DIR}"/*.decoy.faa
        # Study sequences (retired step 1d) in an older release's refs: removed
        # here, offline, with a manifest of what was removed/recovered/lost.
        rm -f "${REFS_DIR}/extra_positives_manifest.tsv"
        python pipeline/strip_ref_tag.py --refs "${REFS_DIR}" --tag _study --threads "${THREADS}"
        # HMM decoys added to families.yaml since the source release: fetch
        # only those (the ones it already holds keep their sequences).
        python pipeline/01_fetch_refs.py --families "${FAMILIES_FILE}" --out "${REFS_DIR}" \
            --only-hmm-decoys
    else
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

    # (Step 1d, which merged local study sequences from extra_sequences/, was
    # retired: those files were not in git, so a rebuild could not reproduce
    # them. Refs from older releases are cleaned in the REUSE_REFS_FROM block.)

    echo ""
    echo "--- 1e. Merge extra positives from NCBI RefSeq (opt-in per family) ---"
    python pipeline/01e_add_refseq_positives.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}"

    echo ""
    echo "--- 1f. Shared random bacterial background (benchmark-only negatives) ---"
    python pipeline/01f_build_random_background.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}" --threads "${THREADS}"

    echo ""
    echo "--- 2. CD-HIT cluster positives + negatives (remove redundancy) ---"
    bash pipeline/02_cluster_cdhit.sh "${REFS_DIR}"
    fi

    echo ""
    echo "--- 3. Train/calib/test split (held-out benchmark set, no leakage; mode=${SPLIT_MODE}) ---"
    python pipeline/03_split_train_test.py --refs "${REFS_DIR}" --families "${FAMILIES_FILE}" \
        --split-mode "${SPLIT_MODE}"

    if [ -n "${REUSE_REFS_FROM:-}" ] && [ "${REUSE_ALIGNMENTS:-1}" = "1" ]; then
        # positive.train is unchanged (step 3 verifies nothing here, but the
        # split is deterministic), so the source release's alignments are
        # the alignments of these sequences. Reusing them keeps the HMMs
        # identical to the source release's: multi-threaded MAFFT is not
        # reproducible run to run (18 of 44 v9 models came out slightly
        # different when realigned from identical input), hmmbuild is.
        echo ""
        echo "--- 4. Reusing alignments from ${REUSE_REFS_FROM} (REUSE_ALIGNMENTS=0 to realign) ---"
        # A family whose train set differs from the source release's (e.g. its
        # study sequences were stripped, which also changes the split) cannot
        # reuse its alignment: only those families are realigned.
        CHANGED=""
        for f in "${REFS_DIR}"/*.positive.train.faa; do
            cmp -s "${f}" "${REUSE_REFS_FROM}/refs/$(basename "${f}")" || \
                CHANGED="${CHANGED} $(basename "${f}" .positive.train.faa)"
        done
        cp "${REUSE_REFS_FROM}"/alignments/* "${ALN_DIR}/"
        if [ -n "${CHANGED}" ]; then
            echo "Train set changed, realigning:${CHANGED}"
            ALIGN_ONLY="${CHANGED# }" THREADS="${THREADS}" \
                bash pipeline/04_align_trim.sh "${REFS_DIR}" "${ALN_DIR}" "${FAMILIES_FILE}"
        fi
    else
    echo ""
    echo "--- 4. Align + trim (TRAIN positives only) ---"
    THREADS="${THREADS}" bash pipeline/04_align_trim.sh "${REFS_DIR}" "${ALN_DIR}" "${FAMILIES_FILE}"
    fi

    echo ""
    echo "--- 5. hmmbuild (or fetch curated Pfam model) + score positive-test/negative sets ---"
    THREADS="${THREADS}" bash pipeline/05_build_hmms.sh "${REFS_DIR}" "${ALN_DIR}" "${HMM_DIR}" "${FAMILIES_FILE}"

    echo ""
    echo "--- 5b. Build HMM decoys (confuser-class models, e.g. peptide ABC ATPases vs opuBA) ---"
    python pipeline/05b_build_hmm_decoys.py --refs "${REFS_DIR}" --aln "${ALN_DIR}" \
        --hmms "${HMM_DIR}" --families "${FAMILIES_FILE}" --threads "${THREADS}" \
        ${REUSE_REFS_FROM:+--reuse-alignments}

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
    CALIB_DIR="${RESULTS_DIR}/calib"
    CASCADE_CONFIG="${RELEASE_DIR}/osmo_refdb.profile_cascade.tsv"
    REPORT_LEVELS="${RELEASE_DIR}/osmo_refdb.profile_report_levels.tsv"
    mkdir -p "${CALIB_DIR}"

    # --- Calibration pass: writes the cascade config and report levels ---
    # Both profile-mode choices are made on reads from positive.calib +
    # negative.train, sequences the benchmark never scores (up to v9 the
    # cascade was tuned on the benchmark reads themselves): which families
    # get an HMM second opinion and at what short-read threshold
    # (osmo_refdb.profile_cascade.tsv), and at which level each family or
    # paralog group is reported (osmo_refdb.profile_report_levels.tsv). The
    # two are decided together per unit: whether a group reaches the bar
    # can depend on whether it gets the cascade.
    echo ""
    echo "--- 9c. Simulate CALIBRATION reads (positive.calib + negative.train) ---"
    python pipeline/09_simulate_reads.py \
        --refs "${REFS_DIR}" --families "${FAMILIES_FILE}" \
        --reads-per-sequence 20 --split calib \
        --out "${CALIB_DIR}/reads"
    python pipeline/09b_compute_read_truth.py --reads "${CALIB_DIR}/reads"

    echo ""
    echo "--- 10c. DIAMOND (cascade scores for every call, none dropped) + HMM on the calibration reads ---"
    CASCADE=score THREADS="${THREADS}" JOBS="${JOBS:-1}" bash pipeline/10_run_benchmark.sh \
        "${CALIB_DIR}/reads" "${CALIB_DIR}" \
        "${RELEASE_DIR}" "${HMM_DIR}/osmo_refdb.hmm"

    echo ""
    echo "--- 11c. Calibration metrics -> report levels + profile cascade config ---"
    python pipeline/11_compute_metrics.py \
        --results "${CALIB_DIR}" --out "${CALIB_DIR}/metrics" --families "${FAMILIES_FILE}" \
        --cascade-out "${CASCADE_CONFIG}" --report-levels-out "${REPORT_LEVELS}"

    # --- Benchmark pass: held-out test reads, cascade on (the default) ---
    echo ""
    echo "--- 9. Simulate reads from held-out test positives, hard negatives + random background ---"
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
    echo "--- 10. osmotool profile (DIAMOND + default cascade) and HMM on the test reads ---"
    CASCADE=on THREADS="${THREADS}" JOBS="${JOBS:-1}" bash pipeline/10_run_benchmark.sh \
        "${RESULTS_DIR}/reads" "${RESULTS_DIR}" \
        "${RELEASE_DIR}" "${HMM_DIR}/osmo_refdb.hmm"

    echo ""
    echo "--- 11. Precision/recall/F1 for DIAMOND, DIAMOND+cascade and HMM + ROC/PR curves ---"
    python pipeline/11_compute_metrics.py \
        --results "${RESULTS_DIR}" --out "${RESULTS_DIR}/metrics" --families "${FAMILIES_FILE}" \
        --cascade-config "${CASCADE_CONFIG}" --report-levels "${REPORT_LEVELS}"

    echo ""
    echo "--- 11b. Annotate-mode benchmark on held-out full-length proteins ---"
    python pipeline/11b_benchmark_annotate.py \
        --refs "${REFS_DIR}" --release "${RELEASE_DIR}" --families "${FAMILIES_FILE}" \
        --out "${RESULTS_DIR}/annotate" --threads "${THREADS}"

    echo ""
    echo "=== Benchmark complete. See ${RESULTS_DIR}/metrics for summary.tsv, "
    echo "    best_threshold_summary.tsv, and ROC/PR plots. ==="
    echo "    Profile-mode DIAMOND+HMM cascade config (from calibration reads): ${CASCADE_CONFIG}"
    echo "    Profile report levels (from calibration reads): ${REPORT_LEVELS}"
    echo "    Per reporting unit on the test reads: ${RESULTS_DIR}/metrics/summary_by_report_unit.tsv"
    echo "    Annotate-mode benchmark: ${RESULTS_DIR}/annotate/summary.tsv"
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
    markers)
        do_markers
        ;;
    nacl_model)
        do_nacl_model
        ;;
    catalog)
        do_catalog
        ;;
    *)
        usage
        ;;
esac
