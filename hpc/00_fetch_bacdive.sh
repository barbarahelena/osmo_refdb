#!/bin/bash
# BacDive strains and their genomes, into <DATA_DIR> (login node, internet; no credentials).
#   1. BacDive harvest (sweep, ~1,820 API calls, about 1.5 h measured): batches already in
#      <DATA_DIR>/bacdive_raw are skipped, so putting the shipped snapshot there first makes
#      this reproduce the published table instead of today's BacDive.
#   2. phenotype table, 3. genomes from NCBI (about 9.4 GB, about 50 min at 8 parallel measured)
#   4. the within-domain salt-in control's strains and genomes (about 0.7 GB)
# Everything is resumable. Run inside the image or any shell with curl, awk and python3+yaml.
#   bash hpc/00_fetch_bacdive.sh data/phenotype_validation
set -euo pipefail
PV="${1:?data dir}"
SCRIPTS="$(cd "$(dirname "$0")/.." && pwd)/pipeline/phenotype_validation"
mkdir -p "${PV}"
bash "${SCRIPTS}/harvest_bacdive.sh" "${PV}" sweep
python3 "${SCRIPTS}/build_table.py" --dir "${PV}"
PAR="${PAR:-8}" bash "${SCRIPTS}/fetch_genomes.sh" "${PV}"
python3 "${SCRIPTS}/build_strategy_control.py" --dir "${PV}"
PAR="${PAR_CONTROL:-4}" bash "${SCRIPTS}/fetch_genomes.sh" "${PV}/strategy_control"
echo "done: $(ls "${PV}/genomes" | wc -l) genomes, $(ls "${PV}/strategy_control/genomes" 2>/dev/null | wc -l) control genomes"
