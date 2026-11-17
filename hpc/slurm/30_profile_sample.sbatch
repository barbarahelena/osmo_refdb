#!/bin/bash
# osmotool profile on one sample, with the catalog's full index if
# FULL_INDEX_DIR exists. Submit per sample:
#   sbatch --export=ALL,R1=...,R2=...,OUT=... hpc/slurm/30_profile_sample.sbatch
# Samples are at most about 20M reads here. Resource requests are placeholders.
#SBATCH --job-name=osmo_profile
#SBATCH --cpus-per-task=8         # ADJUST
#SBATCH --mem=32G                 # ADJUST: the index is loaded in memory
#SBATCH --time=04:00:00           # ADJUST
#SBATCH --output=osmo_profile_%j.log
set -euo pipefail
: "${R1:?set R1}" "${R2:?set R2}" "${OUT:?set OUT (output prefix)}"
cd "${REPO}"
EXTRA=()
[ -f "${FULL_INDEX_DIR}/catalog_index.json" ] && EXTRA+=(--catalog_index "${FULL_INDEX_DIR}")
apptainer exec --pwd /work "${SIF}" osmotool profile "releases/${RELEASE}" \
    -1 "${R1}" -2 "${R2}" --out_prefix "${OUT}" --threads "${SLURM_CPUS_PER_TASK}" \
    --catalog "${CATALOG_NAME}" "${EXTRA[@]+"${EXTRA[@]}"}" --keep_aln
echo "done; sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem"
