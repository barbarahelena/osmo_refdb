#!/bin/bash
# Annotate the gene catalog with the release's HMMs, package it into the
# release (small hits-only index, catalogs/<name>/), then build the full
# competitive index separately. One-off per release (annotation) / catalogue
# (full index). Resource requests are placeholders: run once on the 40 MB
# sample first (hpc/README.md), then set them from sacct's numbers.
#SBATCH --job-name=osmo_catalog_build
#SBATCH --cpus-per-task=16        # ADJUST
#SBATCH --mem=64G                 # ADJUST: bwa-mem3 index build peak memory is the unknown
#SBATCH --time=24:00:00           # ADJUST
#SBATCH --output=osmo_catalog_build_%j.log
set -euo pipefail
cd "${REPO}"
echo "node: $(hostname); AVX2: $(grep -c avx2 /proc/cpuinfo | head -1) cores with avx2"
CATALOG_INDEX=hits THREADS="${SLURM_CPUS_PER_TASK}" CATALOG_FORCE="${CATALOG_FORCE:-0}" \
    apptainer exec --pwd /work "${SIF}" bash run_pipeline.sh catalog "${RELEASE}"
# the full index depends only on the gene file and the aligner, not on the release
apptainer exec --pwd /work "${SIF}" osmotool index-catalog \
    --genes "${CATALOG_GENES}" -o "${FULL_INDEX_DIR}"
echo "done; sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem"
