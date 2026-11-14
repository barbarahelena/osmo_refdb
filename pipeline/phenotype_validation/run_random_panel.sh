#!/usr/bin/env bash
# run_random_panel.sh <outdir> <pool.hmm> [n_parallel]
# Run inside the pinned image (needs prodigal + hmmsearch):
#   docker run --rm -v "$PWD":/work -w /work --entrypoint bash \
#     osmo_refdb:latest pipeline/phenotype_validation/run_random_panel.sh \
#     releases/v8/phenotype_validation_full \
#     releases/v8/phenotype_validation_full/random_panel/pool.hmm 8
#
# Per genome in <outdir>/genomes/, writes <outdir>/random_panel/strains/<bid>/:
#   <bid>.genome_stats.tsv  genome_bp, gc, n_contigs, total_proteins
#                           (the size/GC baselines)
#   <bid>.pool_counts.tsv   proteins hit per pool family (hmmsearch --cut_ga)
#   <bid>.faa.gz            Prodigal proteins, kept gzipped so a future
#                           baseline never needs Prodigal rerun
#
# Prodigal runs with -p meta, exactly as osmotool annotate does, so the
# protein set searched here is the one the osmo panel was searched against.
# One thread per job, parallel across genomes (see run_annotate.sh).
set -uo pipefail
export PATH=/opt/conda/envs/osmo_refdb/bin:$PATH
PV=${1:?usage: run_random_panel.sh <outdir> <pool.hmm> [par]}
POOL=${2:?missing pool.hmm}
PAR=${3:-8}
OUT="$PV/random_panel/strains"
mkdir -p "$OUT"

one() {
  g="$1"; OUT="$2"; POOL="$3"
  b=$(basename "$g"); bid=${b%%_*}
  d="$OUT/$bid"
  [ -s "$d/$bid.pool_counts.tsv" ] && { echo "$bid cached"; return 0; }
  mkdir -p "$d"
  t0=$SECONDS
  tmp="$d/$bid.input.fna"
  case "$g" in *.gz) gunzip -c "$g" > "$tmp" ;; *) cp "$g" "$tmp" ;; esac || {
    echo "$bid rc=input_failed"; rm -f "$tmp"; return 0; }

  # genome size counts A/C/G/T only (N runs and IUPAC codes excluded)
  stats=$(awk 'BEGIN{n=0; gc=0; at=0}
       /^>/ {n++; next}
       { s=toupper($0); gc += gsub(/[GC]/, "", s); at += gsub(/[AT]/, "", s) }
       END{ printf "%d\t%.5f\t%d", gc+at, (gc+at)?gc/(gc+at):0, n }' "$tmp")

  faa="$d/$bid.faa"
  if ! prodigal -i "$tmp" -a "$faa" -o /dev/null -p meta -q 2>"$d/$bid.log"; then
    echo "$bid rc=prodigal_failed"; rm -f "$tmp" "$faa"; return 0
  fi
  rm -f "$tmp"
  if ! hmmsearch --cut_ga --cpu 1 --noali --tblout "$d/$bid.tbl" -o /dev/null \
       "$POOL" "$faa" 2>>"$d/$bid.log"; then
    echo "$bid rc=hmmsearch_failed"; rm -f "$faa"; return 0
  fi
  # proteins hit per family; a protein can hit several families, counted once each
  awk '!/^#/ {k=$3 "\t" $1; if (!(k in seen)) {seen[k]=1; c[$3]++}}
       END{print "family\tn_proteins"; for (f in c) print f "\t" c[f]}' \
    "$d/$bid.tbl" > "$d/$bid.pool_counts.tmp"
  nprot=$(grep -c '^>' "$faa")
  printf 'genome_bp\tgc\tn_contigs\ttotal_proteins\n%s\t%s\n' "$stats" "$nprot" \
    > "$d/$bid.genome_stats.tsv"
  gzip -f "$faa"; rm -f "$d/$bid.tbl"
  # pool_counts.tsv is the "done" marker checked above, so it lands last
  mv "$d/$bid.pool_counts.tmp" "$d/$bid.pool_counts.tsv"
  echo "$bid rc=0 wall=$((SECONDS-t0))s proteins=$nprot"
}
export -f one

T0=$SECONDS
find "$PV/genomes/" -maxdepth 1 \( -name '*.fna' -o -name '*.fna.gz' \) | sort > "$OUT/.inputs"
xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {} "$OUT" "$POOL" < "$OUT/.inputs"
N=$(wc -l < "$OUT/.inputs" | tr -d ' '); W=$((SECONDS-T0))
echo "=== ${PAR}x1: ${W}s wall for $N genomes ==="
