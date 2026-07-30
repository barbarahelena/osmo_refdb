#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 08_build_diamond_db.sh
# Build a DIAMOND database from each family's TRAIN positives (never the
# TEST split, to avoid leakage into the benchmark), plus:
#   * any decoy references from 08a_build_decoy_refs.py (families.yaml:
#     decoy_from_negatives) -- confusable-paralog sequences that must be
#     searchable so they can win select_best_hits' contest away from a
#     mislabeled call, but are excluded from all reported output (see
#     08c_write_scope_manifest.py).
#   * any fused-ORF references from 08d_build_fusion_refs.py (families.yaml:
#     fusion_partner) -- e.g. mrpA_mrpD_fused -- real, reportable detection
#     targets, unlike decoys, so NOT excluded via 08c.
#
# Deliberately EXCLUDES "_study"/"_refseq"-tagged sequences (01d/01e's
# Bakta/RefSeq merges) from the DIAMOND reference specifically, even though
# they remain in positive.train.faa for 04_align_trim.sh's alignment (which
# already ran on the full file before this step, feeding HMM). Confirmed via
# a clean A/B benchmark across 4 families (mrpG, otsA, mscS, otsB -- see
# docs/CHANGELOG.md): the merged pool's added divergent sequences
# consistently helped HMM's profile-based matching (F1 improved in all 4)
# but mostly hurt DIAMOND's identity-based best-hit search (F1 regressed in
# 3/4, up to -0.123) -- DIAMOND does better with a tighter, UniProt-curated
# reference even though HMM benefits from the broader one. Negatives are
# untouched: there's no RefSeq/study-equivalent negative source (01d/01e
# only ever touch positives), so both methods still share the same
# UniProt-only negative pool.
#
# Run 08a and 08d before this script.
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
: > "${COMBINED_FAA}"
for TRAIN_FAA in "${REFS_DIR}"/*.positive.train.faa; do
    [ -s "${TRAIN_FAA}" ] || continue
    # Drop whole records (header + sequence lines) whose tag is _study or
    # _refseq -- `keep` is set once per header and holds for every sequence
    # line until the next one.
    awk '
        /^>/ { keep = ($0 !~ /_study\||_refseq\|/) }
        keep
    ' "${TRAIN_FAA}" >> "${COMBINED_FAA}"
done
if compgen -G "${REFS_DIR}/*.decoy.faa" > /dev/null; then
    cat "${REFS_DIR}"/*.decoy.faa >> "${COMBINED_FAA}"
fi
if compgen -G "${REFS_DIR}/*_fused.faa" > /dev/null; then
    cat "${REFS_DIR}"/*_fused.faa >> "${COMBINED_FAA}"
fi

diamond makedb --in "${COMBINED_FAA}" --db "${RELEASE_DIR}/${RELEASE_NAME}" --threads "${THREADS}"

echo "Built DIAMOND db: ${RELEASE_DIR}/${RELEASE_NAME}.dmnd (excludes _study/_refseq-tagged sequences)"
