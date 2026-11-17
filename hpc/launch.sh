#!/bin/bash
# One script to launch the whole canonical run on the cluster (SLURM + Apptainer).
#
#   cp hpc/slurm/env.sh.example hpc/slurm/env.sh      # edit once
#   source hpc/slurm/env.sh
#   export RELEASE=v10-canonical REUSE_REFS_FROM=releases/v10   # the snapshot to build from
#
#   bash hpc/launch.sh fetch     # login node, internet: GMGC human gut + BacDive strains/genomes
#   bash hpc/launch.sh build     # submit: full index | BacDive proteins | release build ->
#                                #         catalog annotation -> screen          (stops at the gate)
#   ... read the screen (releases/$RELEASE/catalogs/$CATALOG_NAME/screen.tsv) ...
#   bash hpc/launch.sh analyze   # after the gate: BacDive annotation -> analysis + NaCl model
#   bash hpc/launch.sh status
#
# Needs in place before `build` (not downloaded by this script; rsync them, see hpc/README.md):
#   - the image ($SIF), built once from hpc/apptainer/osmo_refdb_hpc.def
#   - the release snapshot at $REUSE_REFS_FROM (refs/, alignments/, markers/, ...)
#   - $DATA_DIR/random_panel/pool.hmm (+ pool_families.tsv, Pfam.version) and, if you want the
#     published BacDive table rather than today's, $DATA_DIR/bacdive_raw/
set -euo pipefail

need() { for v in "$@"; do [ -n "${!v:-}" ] || { echo "ERROR: $v is not set (source hpc/slurm/env.sh)" >&2; exit 1; }; done; }
ok_file() { [ -s "$1" ] || { echo "ERROR: missing $1  ($2)" >&2; exit 1; }; }
REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"; export REPO
cd "${REPO}"
HPC="${REPO}/hpc"

stage="${1:-}"
case "${stage}" in
fetch)
    need CATALOG_GENES DATA_DIR
    echo "### fetch: GMGC human gut -> $(dirname "${CATALOG_GENES}")"
    if [ ! -s "${CATALOG_GENES}" ]; then
        bash "${HPC}/00_fetch_gmgc.sh" "$(dirname "${CATALOG_GENES}")"
    else
        echo "already there: ${CATALOG_GENES}"
    fi
    echo "### fetch: BacDive strains and genomes -> ${DATA_DIR}"
    bash "${HPC}/00_fetch_bacdive.sh" "${DATA_DIR}"
    echo "### fetch done. Pfam is only needed to re-derive the random-panel pool: bash hpc/00_fetch_pfam.sh ${DATA_DIR}"
    ;;

build)
    need SIF RELEASE REUSE_REFS_FROM CATALOG_NAME CATALOG_GENES FULL_INDEX_DIR DATA_DIR NCHUNKS
    ok_file "${SIF}" "the Apptainer image"
    ok_file "${REUSE_REFS_FROM}/refs/manifest.tsv" "the release snapshot to build from"
    ok_file "${CATALOG_GENES}" "run: bash hpc/launch.sh fetch"
    ok_file "${DATA_DIR}/random_panel/pool.hmm" "the random-panel pool snapshot"
    [ -d "${DATA_DIR}/genomes" ] || { echo "ERROR: no ${DATA_DIR}/genomes (run: bash hpc/launch.sh fetch)" >&2; exit 1; }
    [ ! -e "releases/${RELEASE}" ] || { echo "ERROR: releases/${RELEASE} exists; pick a new RELEASE name" >&2; exit 1; }
    LAST=$((NCHUNKS - 1))
    IDX=$(sbatch --parsable "${HPC}/slurm/11_index_catalog.sbatch")
    PROT=$(sbatch --parsable --array=0-"${LAST}" "${HPC}/slurm/40_bacdive_proteins.sbatch")
    BUILD=$(sbatch --parsable "${HPC}/slurm/05_pipeline.sbatch")
    CAT=$(sbatch --parsable --dependency=afterok:"${BUILD}" "${HPC}/slurm/10_build_catalog.sbatch")
    SCR=$(sbatch --parsable --dependency=afterok:"${CAT}" "${HPC}/slurm/15_catalog_screen.sbatch")
    cat <<EOF
submitted (independent, start at once):  full catalog index ${IDX} | BacDive proteins array ${PROT}
submitted (chained):                     release build ${BUILD} -> catalog annotation ${CAT} -> screen ${SCR}

THE GATE: when the screen job ${SCR} ends, read
  releases/${RELEASE}/catalogs/${CATALOG_NAME}/screen.tsv   (and osmo_catalog_screen_${SCR}.log)
A flagged family means: look at those genes. A real confuser class means a decoy and a rebuild
(then rerun 'build' with a new RELEASE); a false alarm means carry on with:
  bash hpc/launch.sh analyze
EOF
    ;;

analyze)
    need SIF RELEASE DATA_DIR NCHUNKS
    ok_file "releases/${RELEASE}/hmms/osmo_refdb.hmm" "the release must be built (launch.sh build)"
    ok_file "${DATA_DIR}/random_panel/pool.hmm" "snapshot"
    N_GEN=$(ls "${DATA_DIR}/genomes" | grep -c '\.fna' || true)
    N_PROT=$(ls "${DATA_DIR}"/random_panel/strains/*/*.pool_counts.tsv 2>/dev/null | wc -l)
    if [ "${N_PROT}" -lt "${N_GEN}" ]; then
        echo "ERROR: the BacDive proteins job has finished ${N_PROT} of ${N_GEN} genomes (squeue -u \$USER); wait for it" >&2
        exit 1
    fi
    bash "${HPC}/prepare_bacdive_release.sh" "${RELEASE}" "${DATA_DIR}"
    LAST=$((NCHUNKS - 1))
    ANN=$(sbatch --parsable --array=0-"${LAST}" "${HPC}/slurm/41_bacdive_annotate.sbatch")
    FIN=$(sbatch --parsable --dependency=afterok:"${ANN}" "${HPC}/slurm/42_bacdive_analyze_nacl.sbatch")
    echo "submitted: BacDive annotation array ${ANN} -> analysis + NaCl model ${FIN}"
    ;;

status)
    squeue -u "${USER}" -o "%.10i %.28j %.8T %.10M %.6D %R" 2>/dev/null || echo "squeue not available here"
    ;;

*)
    sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
    ;;
esac
