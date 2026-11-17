#!/usr/bin/env bash
# run_annotate.sh <outdir> <release> [n_parallel] [threads_each]
# Run inside the pinned image:
#   docker run --rm -v "$PWD":/work -w /work --entrypoint bash \
#     osmo_refdb:latest bin/phenotype_validation/run_annotate.sh \
#     releases/v8/phenotype_validation releases/v8 10 1
#
# One output dir per strain under <outdir>/annotate/. Parallel by strain with
# ONE thread each: hmmscan got slower as threads rose (13s at --cpu 1 vs 16s
# at --cpu 8 on a 5330-protein genome, osmotool 2826ebc), so job-level
# parallelism is the only kind that helped. osmotool v0.7.0 uses hmmsearch
# instead; Prodigal is single-threaded and sets the per-genome floor.
#
# Genomes may be plain (.fna, the pilot) or gzipped (.fna.gz, fetch_genomes.sh
# now). osmotool can't read gzip, so a gzipped genome is unpacked into the
# strain's output dir for the length of its job and deleted afterwards.
#
# Prodigal's protein FASTA is ~99% of annotate's output (111 of 112 MB on the
# pilot) and nothing downstream reads it, so it is dropped by default. Set
# KEEP_PROTEINS=1 to keep it (adds ~17 GB on the full sweep).
set -uo pipefail
export PATH=/opt/conda/envs/osmo_refdb/bin:$PATH
PV=${1:?usage: run_annotate.sh <outdir> <release> [par] [threads]}
REL=${2:?missing release dir}
PAR=${3:-10}; THR=${4:-1}
export KEEP_PROTEINS=${KEEP_PROTEINS:-0}
mkdir -p "$PV/annotate"

one() {
  g="$1"; THR="$2"; PV="$3"; REL="$4"
  b=$(basename "$g"); bid=${b%%_*}
  d="$PV/annotate/$bid"
  [ -s "$d/$bid.gene_counts.tsv" ] && { echo "$bid cached"; return 0; }
  mkdir -p "$d"
  t0=$SECONDS
  input="$g"
  case "$g" in
    *.gz) input="$d/$bid.input.fna"
          gunzip -c "$g" > "$input" || { echo "$bid rc=gunzip_failed"; rm -f "$input"; return 0; } ;;
  esac
  keep=(); [ "$KEEP_PROTEINS" = 1 ] && keep=(--keep_proteins)
  osmotool annotate "$REL" "$input" --out_prefix "$d/$bid" \
    --threads "$THR" "${keep[@]+"${keep[@]}"}" >"$d/$bid.log" 2>&1
  rc=$?
  [ "$input" != "$g" ] && rm -f "$input"
  echo "$bid rc=$rc wall=$((SECONDS-t0))s"
  [ $rc -ne 0 ] && tail -5 "$d/$bid.log"
  return 0
}
export -f one

T0=$SECONDS
# .fna.gz and the pilot's plain .fna; skip .part files from an interrupted fetch
find "$PV/genomes/" -maxdepth 1 \( -name '*.fna' -o -name '*.fna.gz' \) | sort > "$PV/annotate/.inputs"
xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {} "$THR" "$PV" "$REL" < "$PV/annotate/.inputs"
N=$(wc -l < "$PV/annotate/.inputs" | tr -d ' '); W=$((SECONDS-T0))
echo "=== ${PAR}x${THR}: ${W}s wall for $N genomes = $(python3 -c "print(f'{$W/$N:.2f}')")s/genome ==="
