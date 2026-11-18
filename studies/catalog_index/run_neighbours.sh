#!/bin/bash
# The middle option, after run.sh's jobs: the hits index plus the genes close to a hit gene, profiled
# on the same reads, then the scoring of all four indexes again.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "${HERE}/config.sh"
[ -f "${STUDY_DIR}/catalogs/all/catalog.annotation.tsv" ] || { echo "run run.sh first" >&2; exit 1; }
sub() {   # sub <name> <cpus> <mem> <time> <dependency or ''> <step args...>
    local name=$1 cpus=$2 mem=$3 time=$4 dep=$5; shift 5
    sbatch --parsable -p "${PARTITION}" -J "catidx_${name}" -c "${cpus}" --mem "${mem}" -t "${time}" \
        -o "${STUDY_DIR}/logs/${name}.out" ${dep:+--dependency=afterok:${dep}} \
        --wrap "bash ${HERE}/steps.sh $*"
}
nb=$(sub neighbours "${THREADS}" 32G 3:00:00 '' neighbours)
b_nb=$(sub build_neighbours "${THREADS}" 32G 3:00:00 "${nb}" build neighbours)
profiles=()
for div in ${DIVERGENCES}; do
    profiles+=("$(sub "profile_${div}_neighbours" "${THREADS}" 32G 6:00:00 "${b_nb}" profile "${div}" neighbours)")
done
score=$(sub score "1" 48G 2:00:00 "$(IFS=:; echo "${profiles[*]}")" score)
echo "submitted ${nb} ${b_nb} ${profiles[*]} ${score}"
