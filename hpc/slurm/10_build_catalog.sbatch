#!/bin/bash
# Annotate the gene catalog with the release's HMMs and package it into the
# release (annotation + the small hits-only index, catalogs/<name>/). Needs the
# release built first (05_pipeline.sbatch). The full competitive index does not
# depend on the release: build it at the start, in parallel (11_index_catalog.sbatch).
# Measured per stage on 159,560 real gut genes (133 Mbp), then scaled linearly by the
# 75-86x size of the full file (11.9-13.8 M genes): translation (single-threaded Python)
# 8.0 s -> about 10-12 min; hmmsearch 17.2 s at 4 threads -> about 21-25 min at 4 threads;
# Python memory 203 bytes per gene -> 2.4-2.8 GB; the small hits-only index (about 5% of
# genes) is a few GB at most. Requests: about 2x memory and 4x time.
#SBATCH --job-name=osmo_catalog_annotate
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=osmo_catalog_annotate_%j.log
set -euo pipefail
cd "${REPO}"
echo "node: $(hostname); AVX2: $(grep -c avx2 /proc/cpuinfo | head -1) cores with avx2"
CATALOG_INDEX=hits THREADS="${SLURM_CPUS_PER_TASK}" CATALOG_FORCE="${CATALOG_FORCE:-0}" \
    apptainer exec --pwd /work "${SIF}" bash run_pipeline.sh catalog "${RELEASE}"
echo "done; sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem"
