#!/bin/bash
#SBATCH --job-name=mixstudy_full
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=work/%x-%j.log
# The full index of every catalog gene (built once into FULL_INDEX, reused if there), then
# the same reads mapped to it and scored: work/gmgc_full/, and work/index_comparison.md.
# Submitted by submit.sh after study.sbatch, with --chdir at the study.
#
# Memory and time are NOT measured at this size. 240 GB and 12 h follow the pipeline's
# former OSMOTOOL_INDEXCATALOG request, itself scaled from bwa-mem3 on 133 Mbp of the
# catalog (1.37 GB peak, 33 s). This job records what the index really took:
# work/gmgc_full/index.time (GNU time) and the sacct line at the end of the log.
set -euo pipefail
source hpc/env.sh
export CATALOG_INDEX=$FULL_INDEX CATALOG_TAG=gmgc_full
bash run_study.sh index map score verify
$PYTHON_RUN -I "$STUDY/compare_indexes.py" --work "$STUDY/work" --small gmgc_gut --full gmgc_full
du -sh "$FULL_INDEX"
grep -E "Maximum resident|Elapsed" work/gmgc_full/index.time 2>/dev/null || true
sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,MaxRSS,ReqMem,NCPUS 2>/dev/null || true
