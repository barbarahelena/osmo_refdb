#!/bin/bash
# Submit the catalog-index study to SLURM, each step after the ones it needs (00_setup.sh first).
# Results: ${STUDY_DIR}/results/ (summary.tsv, by_family.tsv, verification.tsv, resources.tsv).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "${HERE}/config.sh"
[ -f "${CATALOG_GENES}" ] || { echo "no catalog at ${CATALOG_GENES} (set CATALOG_GENES)" >&2; exit 1; }
[ -x "${ENV_DIR}/bin/osmotool" ] || { echo "run 00_setup.sh first" >&2; exit 1; }
mkdir -p "${STUDY_DIR}/logs"
sub() {   # sub <name> <cpus> <mem> <time> <dependency or ''> <step args...>
    local name=$1 cpus=$2 mem=$3 time=$4 dep=$5; shift 5
    sbatch --parsable -p "${PARTITION}" -J "catidx_${name}" -c "${cpus}" --mem "${mem}" -t "${time}" \
        -o "${STUDY_DIR}/logs/${name}.out" ${dep:+--dependency=afterok:${dep}} \
        --wrap "bash ${HERE}/steps.sh $*"
}

slice=$(sub slice 1 8G 2:00:00 '' slice)
# the full index is the one whose memory is in question: room to measure it
b_all=$(sub build_all "${THREADS}" 240G 6:00:00 "${slice}" build all)
b_hits=$(sub build_hits "${THREADS}" 32G 3:00:00 "${slice}" build hits)
profiles=()
for div in ${DIVERGENCES}; do
    sim=$(sub "simulate_${div}" 1 24G 3:00:00 "${slice}" simulate "${div}")
    profiles+=("$(sub "profile_${div}_none" "${THREADS}" 32G 6:00:00 "${sim}" profile "${div}" none)")
    profiles+=("$(sub "profile_${div}_hits" "${THREADS}" 32G 6:00:00 "${sim}:${b_hits}" profile "${div}" hits)")
    profiles+=("$(sub "profile_${div}_all" "${THREADS}" 160G 6:00:00 "${sim}:${b_all}" profile "${div}" all)")
done
score=$(sub score 1 48G 2:00:00 "$(IFS=:; echo "${profiles[*]}")" score)
echo "submitted; the last job (score) is ${score}. Logs: ${STUDY_DIR}/logs/"
