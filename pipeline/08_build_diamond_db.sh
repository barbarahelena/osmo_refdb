#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 08_build_diamond_db.sh
# Build a DIAMOND database from each family's TRAIN positives (never the
# TEST split, to avoid leakage into the benchmark).
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

diamond makedb --in "${COMBINED_FAA}" --db "${RELEASE_DIR}/${RELEASE_NAME}" --threads "${THREADS}"

echo "Built DIAMOND db: ${RELEASE_DIR}/${RELEASE_NAME}.dmnd"
