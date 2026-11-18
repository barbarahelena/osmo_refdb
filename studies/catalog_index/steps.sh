#!/bin/bash
# One step of the catalog-index study (run.sh submits them): bash steps.sh <step> [args]
#   slice                 the seeded gene slice of the catalog
#   build <all|hits>      build-catalog of the slice with that index
#   simulate <div>        read pairs from the slice at that divergence
#   profile <div> <arm>   osmotool profile on those reads: arm none (--no_catalog), hits or all
#   score                 the comparison (score.py) and the resources used
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "${HERE}/config.sh"
export PATH="${ENV_DIR}/bin:${PATH}"
S="${STUDY_DIR}"
mkdir -p "${S}/logs"
timed() { local name=$1; shift; /usr/bin/time -v -o "${S}/logs/${name}.time" "$@"; }

case "$1" in
slice)
    python3 "${HERE}/slice_genes.py" --genes "${CATALOG_GENES}" --n "${N_GENES}" --seed "${SEED}" \
        --out "${S}/slice.fna.tmp"
    mv "${S}/slice.fna.tmp" "${S}/slice.fna"
    ;;
build)
    arm=$2
    rm -rf "${S}/catalogs/${arm}"
    mkdir -p "${S}/catalogs"
    timed "build_${arm}" osmotool build-catalog "${RELEASE}" --genes "${S}/slice.fna" \
        -o "${S}/catalogs/${arm}" --index "${arm}" -t "${THREADS}"
    du -sb "${S}/catalogs/${arm}" > "${S}/logs/build_${arm}.size"
    ;;
simulate)
    div=$2
    python3 "${HERE}/simulate_reads.py" --genes "${S}/slice.fna" --pairs "${N_PAIRS}" --divergence "${div}" \
        --error-rate "${ERROR_RATE}" --seed "${SEED}" --out "${S}/reads/${div}"
    ;;
profile)
    div=$2 arm=$3
    out="${S}/runs/${div}/${arm}"
    rm -rf "${out}"
    mkdir -p "${out}"
    case "${arm}" in
        none) cat_args=(--no_catalog) ;;
        *)    cat_args=(--catalog "${S}/catalogs/${arm}") ;;
    esac
    # profile's own defaults (what a user gets), the alignments kept for the per-read rebuild
    timed "profile_${div}_${arm}" osmotool profile "${RELEASE}" \
        -1 "${S}/reads/${div}_R1.fastq" -2 "${S}/reads/${div}_R2.fastq" \
        --out_prefix "${out}/reads" "${cat_args[@]}" --keep_aln --no_markers --threads "${THREADS}"
    ;;
score)
    python3 "${HERE}/score.py" --compute-metrics "${HERE}/../../modules/local/benchmark/metrics/templates/compute_metrics.py" \
        --release "${RELEASE}" --annotation "${S}/catalogs/all/catalog.annotation.tsv" \
        --runs "${S}/runs" --reads "${S}/reads" --out "${S}/results"
    {
        printf 'step\twall_clock\tmax_rss_gb\n'
        for f in "${S}"/logs/*.time; do
            n=$(basename "$f" .time)
            wall=$(sed -n 's/.*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$f")
            rss=$(awk -F': ' '/Maximum resident set size/ { printf "%.2f", $2 / 1048576 }' "$f")
            printf '%s\t%s\t%s\n' "$n" "$wall" "$rss"
        done
    } > "${S}/results/resources.tsv"
    for f in "${S}"/logs/build_*.size; do
        printf '%s\t%.2f GB\n' "$(basename "$f" .size)" "$(awk '{ print $1 / 1e9 }' "$f")"
    done > "${S}/results/catalog_sizes.tsv"
    cat "${S}/results/resources.tsv" "${S}/results/catalog_sizes.tsv"
    ;;
*)
    echo "unknown step: $1" >&2; exit 1 ;;
esac
