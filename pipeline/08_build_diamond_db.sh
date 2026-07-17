#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 08_build_diamond_db.sh
# Build a DIAMOND database from each family's TRAIN positives (never the
# TEST split, to avoid leakage into the benchmark), plus any decoy
# references from 08a_build_decoy_refs.py (families.yaml:
# decoy_from_negatives) -- confusable-paralog sequences that must be
# searchable so they can win select_best_hits' contest away from a
# mislabeled call, but are excluded from all reported output (see
# 08c_write_scope_manifest.py). Run 08a before this script.
#
# Usage: bash 08_build_diamond_db.sh [refs_dir] [release_dir] [release_name]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
RELEASE_DIR="${2:-releases/dev}"
RELEASE_NAME="${3:-osmo_refdb}"
THREADS="${THREADS:-4}"

mkdir -p "${RELEASE_DIR}"

COMBINED_FAA="${RELEASE_DIR}/${RELEASE_NAME}.train_refs.faa"
cat "${REFS_DIR}"/*.positive.train.faa > "${COMBINED_FAA}"
if compgen -G "${REFS_DIR}/*.decoy.faa" > /dev/null; then
    cat "${REFS_DIR}"/*.decoy.faa >> "${COMBINED_FAA}"
fi

diamond makedb --in "${COMBINED_FAA}" --db "${RELEASE_DIR}/${RELEASE_NAME}" --threads "${THREADS}"

echo "Built DIAMOND db: ${RELEASE_DIR}/${RELEASE_NAME}.dmnd"
