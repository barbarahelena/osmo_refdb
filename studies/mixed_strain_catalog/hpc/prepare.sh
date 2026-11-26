#!/usr/bin/env bash
# On a login node (needs internet): the two container images, the catalog's gene file
# (if missing) and the study's genomes. Then: bash hpc/submit.sh
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p "$SIF_DIR"

[ -s "$OSMOTOOL_SIF" ] || $APPTAINER pull "$OSMOTOOL_SIF" docker://docker.io/barbarahelena/osmotool:0.11.1
[ -s "$WGSIM_SIF" ] || curl -sSfL -o "$WGSIM_SIF" https://depot.galaxyproject.org/singularity/wgsim:1.0--h577a1d6_10

[ -s "$CATALOG/catalog.json" ] || { echo "no catalog in $CATALOG: set RELEASE in hpc/env.sh" >&2; exit 1; }
if [ ! -s "$CATALOG_GENES" ]; then   # 2.99 GB, one member of a 19.3 GB zip on Zenodo, fetched alone
    $PYTHON_RUN -I "$ROOT/modules/local/catalog/fetch/templates/fetch_zenodo_member.py" \
        --member "gene_catalogues/$(basename "$CATALOG_GENES")" --out "$CATALOG_GENES"
fi

bash "$STUDY/run_study.sh" genomes
echo "prepared: images in $SIF_DIR, genes $CATALOG_GENES, genomes in $STUDY/work/genomes"
