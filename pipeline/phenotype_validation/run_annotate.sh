#!/usr/bin/env bash
# run_annotate.sh <outdir> <release> [n_parallel] [threads_each]
# Run inside the pinned image:
#   docker run --rm -v "$PWD":/work -w /work --entrypoint bash \
#     osmo_refdb:latest pipeline/phenotype_validation/run_annotate.sh \
#     releases/v8/phenotype_validation releases/v8 10 1
#
# One output dir per strain under <outdir>/annotate/. Parallel by strain with
# ONE thread each: hmmscan gets slower as threads rise (13s at --cpu 1 vs 16s
# at --cpu 8 on a 5330-protein genome), so job-level parallelism is the only
# kind that helps. Prodigal is single-threaded and sets the per-genome floor.
set -uo pipefail
export PATH=/opt/conda/envs/osmo_refdb/bin:$PATH
PV=${1:?usage: run_annotate.sh <outdir> <release> [par] [threads]}
REL=${2:?missing release dir}
PAR=${3:-10}; THR=${4:-1}
mkdir -p "$PV/annotate"

one() {
  g="$1"; THR="$2"; PV="$3"; REL="$4"
  b=$(basename "$g" .fna); bid=${b%%_*}
  d="$PV/annotate/$bid"
  [ -s "$d/$bid.gene_counts.tsv" ] && { echo "$bid cached"; return 0; }
  mkdir -p "$d"
  t0=$SECONDS
  osmotool annotate "$REL" "$g" --out_prefix "$d/$bid" \
    --threads "$THR" --keep_proteins >"$d/$bid.log" 2>&1
  rc=$?
  echo "$bid rc=$rc wall=$((SECONDS-t0))s"
  [ $rc -ne 0 ] && tail -5 "$d/$bid.log"
  return 0
}
export -f one

T0=$SECONDS
ls "$PV"/genomes/*.fna | xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {} "$THR" "$PV" "$REL"
N=$(ls "$PV"/genomes/*.fna | wc -l | tr -d ' '); W=$((SECONDS-T0))
echo "=== ${PAR}x${THR}: ${W}s wall for $N genomes = $(python3 -c "print(f'{$W/$N:.2f}')")s/genome ==="
