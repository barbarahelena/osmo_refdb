#!/bin/bash
# Download the human-gut GMGC catalog (2.99 GB) from Zenodo, without the
# other 16 GB of the zip. Run on a login/transfer node with internet.
#   bash hpc/00_fetch_gmgc.sh /path/to/dest_dir            # the full file
#   bash hpc/00_fetch_gmgc.sh /path/to/dest_dir 40000000   # first 40 MB only: a smoke-test sample
set -euo pipefail
DEST="${1:?destination directory}"
MAXB="${2:-}"
MEMBER=gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
mkdir -p "${DEST}"
python3 "$(dirname "$0")/fetch_zenodo_member.py" --member "${MEMBER}" \
    --out "${DEST}/$(basename "${MEMBER}")${MAXB:+.head}" ${MAXB:+--max-bytes "${MAXB}"}
