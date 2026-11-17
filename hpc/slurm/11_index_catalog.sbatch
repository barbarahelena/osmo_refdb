#!/bin/bash
# The full, competitive index of the catalog's genes. It depends only on the gene file and
# the aligner, not on any release, so it starts at once and runs alongside the release
# build. Built with CLASSIC bwa, not bwa-mem3 (Cayman does the same), because of memory.
# Measured on 133 Mbp of the real gut catalog:
#   bwa-mem3 index  33 s, peak 1,371 MB (10.3 bytes/bp), 472 MB on disk
#   classic bwa     92 s, peak   203 MB ( 1.5 bytes/bp), 240 MB on disk, single-threaded
# The full file is about 10-11.5 Gbp (75-86x the sample), so linearly scaled: classic bwa
# about 15-17 GB of memory, about 2 h, 18-21 GB on disk. bwa-mem3's own documentation gives a
# higher per-base cost above about 1 Gbp (72 GB for human), which would put it at 100 GB or
# well above: not used for the full index. (The release's small hits-only index stays bwa-mem3.)
#SBATCH --job-name=osmo_catalog_index
#SBATCH --cpus-per-task=1         # bwa index is single-threaded
#SBATCH --mem=32G                 # 2x the scaled estimate
#SBATCH --time=06:00:00           # about 2.5x the scaled estimate
#SBATCH --output=osmo_catalog_index_%j.log
set -euo pipefail
cd "${REPO}"
apptainer exec --pwd /work "${SIF}" osmotool index-catalog --aligner bwa \
    --genes "${CATALOG_GENES}" -o "${FULL_INDEX_DIR}"
echo "done; sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem"
