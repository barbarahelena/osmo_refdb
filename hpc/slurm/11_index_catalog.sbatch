#!/bin/bash
# The full, competitive bwa-mem3 index of the catalog's genes. It depends only
# on the gene file and the aligner, not on any release, so it can start at once
# and run alongside the release build. Resource requests are placeholders.
#SBATCH --job-name=osmo_catalog_index
#SBATCH --cpus-per-task=16        # ADJUST
#SBATCH --mem=64G                 # ADJUST: bwa-mem3 index build memory at full size is not measured
#SBATCH --time=24:00:00           # ADJUST
#SBATCH --output=osmo_catalog_index_%j.log
set -euo pipefail
cd "${REPO}"
apptainer exec --pwd /work "${SIF}" osmotool index-catalog --genes "${CATALOG_GENES}" -o "${FULL_INDEX_DIR}"
echo "done; sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem"
