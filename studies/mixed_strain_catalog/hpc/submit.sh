#!/usr/bin/env bash
# Submit the two jobs: the shipped index (study.sbatch), then the full index (full.sbatch).
#   PARTITION       partition for both (default genoa, as the pipeline's snellius profile)
#   PARTITION_FULL  partition for the full-index job (default PARTITION); e.g. fat_genoa if
#                   240 GB is not enough, after raising --mem in full.sbatch
set -euo pipefail
STUDY=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$STUDY/work"
P=${PARTITION:-genoa}
j1=$(sbatch --parsable --chdir="$STUDY" -p "$P" "$STUDY/hpc/study.sbatch")
j2=$(sbatch --parsable --chdir="$STUDY" -p "${PARTITION_FULL:-$P}" --dependency=afterok:"$j1" "$STUDY/hpc/full.sbatch")
echo "submitted $j1 (shipped index) and $j2 (full index, after $j1); logs in $STUDY/work/"
