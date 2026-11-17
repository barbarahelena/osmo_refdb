#!/bin/bash
# Submit the canonical chain, with the dependencies between the steps, up to the
# gate. After the gate (human review of the screen), BacDive and the NaCl model
# run on the frozen HMMs.
#
#   source hpc/slurm/env.sh
#   export RELEASE=v10-canonical REUSE_REFS_FROM=releases/v10
#   bash hpc/slurm/submit_all.sh
#
#   index (gene file only) --------------------------------------+
#   build the release (HMMs, DIAMOND, benchmarks) -> annotate catalog -> screen
#                                                                         |
#                           GATE: flagged family? yes -> decoy, rebuild; no -> freeze
set -euo pipefail
: "${REPO:?}" "${SIF:?}" "${RELEASE:?}" "${CATALOG_NAME:?}" "${CATALOG_GENES:?}" "${FULL_INDEX_DIR:?}"
cd "${REPO}"
IDX=$(sbatch --parsable hpc/slurm/11_index_catalog.sbatch)
BUILD=$(sbatch --parsable hpc/slurm/05_pipeline.sbatch)
CAT=$(sbatch --parsable --dependency=afterok:"${BUILD}" hpc/slurm/10_build_catalog.sbatch)
SCR=$(sbatch --parsable --dependency=afterok:"${CAT}" hpc/slurm/15_catalog_screen.sbatch)
echo "index ${IDX} | build ${BUILD} -> catalog ${CAT} -> screen ${SCR}"
echo "When the screen job ends: read its log (and ${REPO}/releases/${RELEASE}/catalogs/${CATALOG_NAME}/screen.tsv)."
