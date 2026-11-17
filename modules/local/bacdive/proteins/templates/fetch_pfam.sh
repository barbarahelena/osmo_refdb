#!/bin/bash
# Pfam-A, pinned to release 39.0 (the one the BacDive random-panel pool was drawn
# from), into <DATA_DIR>/random_panel/. Only needed to re-derive the pool
# (build_random_pool.py); the pool itself (pool.hmm) travels as a snapshot. Run on a
# login node (internet). Note: re-deriving the pool also needs screen_proteomes.faa,
# which no script in this repo creates.
#   bash bin/fetch_pfam.sh data/phenotype_validation
set -euo pipefail
DEST="${1:?destination data dir}/random_panel"
BASE=https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam39.0
mkdir -p "${DEST}"
for f in Pfam-A.hmm.gz Pfam.version.gz; do
    [ -s "${DEST}/${f}" ] || curl -sS -f -L --retry 3 -o "${DEST}/${f}.part" "${BASE}/${f}" && mv "${DEST}/${f}.part" "${DEST}/${f}"
done
[ -s "${DEST}/Pfam.version" ] || gunzip -c "${DEST}/Pfam.version.gz" > "${DEST}/Pfam.version"
cat "${DEST}/Pfam.version"
