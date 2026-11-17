#!/bin/bash
# The Nextflow head job: one small SLURM job that submits and watches all the others, so
# nothing runs on the login node. It lives as long as the workflow, so ask for the longest
# walltime your site allows (the wait for the release build and the downloads dominate).
#
# From the repo root:
#   sbatch run_nextflow.sbatch --release v11
# after reading the gate (first run's output):
#   sbatch run_nextflow.sbatch --release v11 -resume --gate_passed true
#
# At a site with an nf-core/configs profile, add  -profile <institution>  after the script name
# (without one the generic SLURM + Apptainer profile is used).
#
# Nextflow must be on the node: `module load nextflow` or put it on PATH below.
#SBATCH --job-name=osmo_nf_head
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=osmo_nf_%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load nextflow 2>/dev/null || true
export NXF_OPTS="${NXF_OPTS:--Xms256m -Xmx2g}"
export NXF_ANSI_LOG=false
# one trace/timeline/report per invocation: the measured peak memory and time of every task
# (what to tighten the requests in nextflow.config with)
INFO="pipeline_info/$(date +%Y%m%d-%H%M%S)_${SLURM_JOB_ID}"
mkdir -p pipeline_info
nextflow run main.nf -with-trace "${INFO}.trace.tsv" -with-timeline "${INFO}.timeline.html" \
    -with-report "${INFO}.report.html" "$@"
