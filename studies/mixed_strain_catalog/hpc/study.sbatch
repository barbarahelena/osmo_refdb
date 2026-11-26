#!/bin/bash
#SBATCH --job-name=mixstudy_shipped
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=work/%x-%j.log
# The steps that do not depend on the index (annotation, reads, truth, DIAMOND), then the
# release's shipped index: work/gmgc_gut/. Submitted by submit.sh with --chdir at the study.
# Resources: locally the whole study ran in minutes on 8 CPUs in a Docker VM of about 7 GB,
# with osmotool emulated (amd64 on arm64); 32 GB and 2 h leave room.
set -euo pipefail
source hpc/env.sh
bash run_study.sh annotate reads truth diamond map score verify

# the reads and truth the README's results were computed from; if they differ (another
# wgsim build), the shipped-index numbers here differ from the README's, but the
# comparison with the full index is still on the same reads
(cd work && md5sum -c ../expected.md5) || echo "WARNING: reads or truth differ from the README's run"
