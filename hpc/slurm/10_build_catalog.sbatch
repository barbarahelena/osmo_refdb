#!/bin/bash
# Annotate the gene catalog with the release's HMMs and package it into the
# release (annotation + the small hits-only index, catalogs/<name>/). Needs the
# release built first (05_pipeline.sbatch). The full competitive index does not
# depend on the release: build it at the start, in parallel (11_index_catalog.sbatch).
# Resource requests are placeholders: run once on the 40 MB sample first
# (hpc/README.md), then set them from sacct's numbers.
#SBATCH --job-name=osmo_catalog_annotate
#SBATCH --cpus-per-task=16        # ADJUST
#SBATCH --mem=64G                 # ADJUST
#SBATCH --time=24:00:00           # ADJUST
#SBATCH --output=osmo_catalog_annotate_%j.log
set -euo pipefail
cd "${REPO}"
echo "node: $(hostname); AVX2: $(grep -c avx2 /proc/cpuinfo | head -1) cores with avx2"
CATALOG_INDEX=hits THREADS="${SLURM_CPUS_PER_TASK}" CATALOG_FORCE="${CATALOG_FORCE:-0}" \
    apptainer exec --pwd /work "${SIF}" bash run_pipeline.sh catalog "${RELEASE}"
echo "done; sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem"
