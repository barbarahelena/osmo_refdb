#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — run_benchmark.sh
# Run DIAMOND (via osmotool profile, train-only db; with the DIAMOND+HMM
# cascade unless CASCADE=off) and HMM (translate + hmmscan, train-only HMMs) on the simulated reads, for later
# scoring by compute_metrics.py.
#
# Usage: bash run_benchmark.sh [reads_dir] [results_dir] [release_dir] [hmm_db]
#
# release_dir is an unpacked osmo_refdb release directory (e.g. releases/v8),
# NOT the .dmnd file path directly -- `osmotool profile` finds osmo_refdb.dmnd,
# hmms/osmo_refdb.hmm, osmo_refdb.profile_cascade.tsv, and
# osmo_refdb.profile_excluded_families.txt inside it by their fixed names.
# ---------------------------------------------------------------------------

set -euo pipefail

READS_DIR="${1:-results/reads}"
OUT_DIR="${2:-results}"
RELEASE_DIR="${3:-releases/dev}"
HMM_DB="${4:-hmms/osmo_refdb.hmm}"
THREADS="${THREADS:-4}"
# The benchmark's positive test set is deliberately broad/phylogenetically
# diverse (all UniProt hits for the gene, not just close orthologs), so
# many true positives align at only 50-70% identity to their nearest
# training-set reference. osmotool profile's published default (--id 80)
# is tuned for high-confidence calls in typical use, but would silently
# drop most divergent true positives here. Relax to more permissive
# defaults for benchmarking; override via env vars if needed.
DIAMOND_MIN_IDENTITY="${DIAMOND_MIN_IDENTITY:-0.40}"
DIAMOND_MIN_QUERY_COVER="${DIAMOND_MIN_QUERY_COVER:-0.50}"
# Off by default (0.0), same as osmotool's own default -- set this to test
# whether filtering on subject/target coverage reduces the multi-domain
# fusion-protein false positives found while benchmarking proX (see
# --min_subject_cover's help text in osmotool for why this isn't on by
# default: a short read is always much shorter than most full-length
# reference proteins, so a strict threshold here also costs real recall).
DIAMOND_MIN_SUBJECT_COVER="${DIAMOND_MIN_SUBJECT_COVER:-0.0}"
# CASCADE=on (default) runs `osmotool profile` the way users get it: with the
# DIAMOND+HMM cascade, configured by <release>/osmo_refdb.profile_cascade.tsv.
# That file is written by step 11 from the CALIBRATION reads (the pipeline
# runs steps 9-11 on --split calib first), never from the reads scored here.
# --keep_aln retains the raw DIAMOND alignments and the cascade's own
# hmmsearch table, so compute_metrics.py scores DIAMOND alone and
# DIAMOND+cascade from this single run.
# CASCADE=score is the calibration run's mode: every family in the HMM db is
# flagged with a threshold no read can fail, so osmotool keeps every DIAMOND
# call (gene_counts.tsv = DIAMOND-only) but its --keep_aln
# <sample>.cascade.hmmsearch.tblout holds each call's HMM score against the
# family it was called to. compute_metrics.py replays any cascade config
# from that, which is how report level and cascade get decided together on
# the calibration reads. The scores are what production would see: the
# cascade's hmmsearch fixes -Z to the number of models, so neither a score
# nor whether it is reported depends on which reads were searched. osmotool
# finds its cascade config by fixed name in the release dir, so this mode
# points it at a shadow dir: symlinks to the release's files plus the
# flag-every-family config.
# CASCADE=off passes --no_cascade: for osmotool versions older than ab3fb9e,
# which have no cascade default and no such flag.
CASCADE="${CASCADE:-on}"
# HMM samples run JOBS at a time, one hmmscan thread each: hmmscan against
# this small db gets slower, not faster, with more threads per job.
JOBS="${JOBS:-1}"
CASCADE_ARGS=()
PROFILE_DB="${RELEASE_DIR}"
if [ "${CASCADE}" = "off" ]; then
    if osmotool profile --help 2>/dev/null | grep -q -- "--no_cascade"; then
        CASCADE_ARGS+=(--no_cascade)
    fi
elif [ "${CASCADE}" = "score" ]; then
    PROFILE_DB="${OUT_DIR:?}/cascade_score_db"
    # only symlinks and one small TSV live here; rebuilt each run
    rm -f "${PROFILE_DB}"/*
    mkdir -p "${PROFILE_DB}"
    release_abs="$(cd "${RELEASE_DIR}" && pwd)"
    for f in "${release_abs}"/*; do
        case "$(basename "${f}")" in
            osmo_refdb.profile_cascade.tsv|osmo_refdb.profile_report_levels.tsv|results|markers) ;;
            *) ln -s "${f}" "${PROFILE_DB}/$(basename "${f}")" ;;
        esac
    done
    {
        printf 'family\tdiamond_precision\tprecision_metric\tneeds_cascade_check\thmm_short_read_threshold\n'
        awk '$1 == "NAME" {printf "%s\t\tscore_all\tyes\t-1000000000\n", $2}' "${HMM_DB}"
    } > "${PROFILE_DB}/osmo_refdb.profile_cascade.tsv"
elif [ "${CASCADE}" != "on" ]; then
    echo "ERROR: CASCADE must be on, off or score (got '${CASCADE}')" >&2
    exit 1
elif [ ! -s "${RELEASE_DIR}/osmo_refdb.profile_cascade.tsv" ]; then
    echo "ERROR: CASCADE=on but ${RELEASE_DIR}/osmo_refdb.profile_cascade.tsv is missing." >&2
    echo "       Run the calibration benchmark first (BENCHMARK_METRICS_CALIB writes it)." >&2
    exit 1
fi

# osmotool >= 0.8.0 runs a second DIAMOND pass against <release>/markers/ for
# copies-per-genome normalisation. The benchmark reads come from single genes,
# not genomes, so there is nothing to normalise to, and the pass does not
# change which reads are counted: skip it.
if osmotool profile --help 2>/dev/null | grep -q -- "--no_markers"; then
    CASCADE_ARGS+=(--no_markers)
fi

# One sample's HMM pass: R1 and R2 together, OrfM's 6-frame ORFs, hmmscan against the HMM db.
# No --cut_ga: the GA cutoffs were calibrated on full-length proteins, which a short read's ORF
# fragment never reaches; compute_metrics.py sweeps the raw bit scores instead, as for DIAMOND.
translate_and_hmmscan() (   # <R1> <R2> <hmm_db> <out_prefix> <threads>
    set -euo pipefail
    local_prefix="$4"
    mkdir -p "$(dirname "${local_prefix}")"
    zcat -f "$1" "$2" > "${local_prefix}.merged.fastq"
    orfm "${local_prefix}.merged.fastq" > "${local_prefix}.orfs.faa"
    hmmscan --cpu "${5:-4}" --tblout "${local_prefix}.hmmscan.tblout" --noali "$3" "${local_prefix}.orfs.faa" > /dev/null
    rm -f "${local_prefix}.merged.fastq"
    echo "  HMM hits: ${local_prefix}.hmmscan.tblout"
)
export -f translate_and_hmmscan

mkdir -p "${OUT_DIR}/diamond" "${OUT_DIR}/hmm"

echo "=== DIAMOND (osmotool profile, cascade ${CASCADE}, --min_identity ${DIAMOND_MIN_IDENTITY} --min_query_cover ${DIAMOND_MIN_QUERY_COVER} --min_subject_cover ${DIAMOND_MIN_SUBJECT_COVER}) ==="
for R1 in "${READS_DIR}"/*_R1.fastq*; do
    [ -e "${R1}" ] || continue
    SAMPLE=$(basename "${R1}" | sed -E 's/_R1\.fastq(\.gz)?$//')
    R2="${READS_DIR}/${SAMPLE}_R2.fastq"
    [ -f "${R2}.gz" ] && R2="${R2}.gz"
    [ -f "${R2}" ] || { echo "  SKIP ${SAMPLE}: no R2 found"; continue; }

    echo "  sample: ${SAMPLE}"
    osmotool profile \
        "${PROFILE_DB}" \
        -1 "${R1}" -2 "${R2}" \
        --out_prefix "${OUT_DIR}/diamond/${SAMPLE}" \
        --threads "${THREADS}" \
        --min_identity "${DIAMOND_MIN_IDENTITY}" \
        --min_query_cover "${DIAMOND_MIN_QUERY_COVER}" \
        --min_subject_cover "${DIAMOND_MIN_SUBJECT_COVER}" \
        --keep_aln \
        "${CASCADE_ARGS[@]+"${CASCADE_ARGS[@]}"}"
done

if [ "${SKIP_HMM:-0}" = "1" ]; then
    echo "SKIP_HMM=1: not running the HMM pass (reusing ${OUT_DIR}/hmm as it is)"
    exit 0
fi

echo ""
echo "=== HMM (6-frame translate + hmmscan; ${JOBS} sample(s) at a time) ==="
HMM_THREADS="${THREADS}"
[ "${JOBS}" -gt 1 ] && HMM_THREADS=1
for R1 in "${READS_DIR}"/*_R1.fastq*; do
    [ -e "${R1}" ] || continue
    SAMPLE=$(basename "${R1}" | sed -E 's/_R1\.fastq(\.gz)?$//')
    R2="${READS_DIR}/${SAMPLE}_R2.fastq"
    [ -f "${R2}.gz" ] && R2="${R2}.gz"
    [ -f "${R2}" ] || { echo "  SKIP ${SAMPLE}: no R2 found" >&2; continue; }
    printf '%s\t%s\t%s\n' "${R1}" "${R2}" "${OUT_DIR}/hmm/${SAMPLE}"
done | xargs -P "${JOBS}" -L 1 bash -c \
    'translate_and_hmmscan "$3" "$4" "$1" "$5" "$2"' \
    _ "${HMM_DB}" "${HMM_THREADS}"

echo ""
echo "Done. Now run: python compute_metrics.py --results ${OUT_DIR} --out ${OUT_DIR}/metrics"
