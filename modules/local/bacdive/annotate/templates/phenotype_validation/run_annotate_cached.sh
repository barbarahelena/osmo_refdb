#!/usr/bin/env bash
# run_annotate_cached.sh <outdir> <release> [n_parallel]
# Like run_annotate.sh, but passes the Prodigal proteins cached by
# run_random_panel.sh (random_panel/strains/<bid>/<bid>.faa.gz, Prodigal -p meta)
# via --proteins, so Prodigal isn't rerun. Checked on 3 genomes against
# releases/v8: gene_counts, systems and gene_coordinates byte-identical to the
# original v8 run_annotate.sh output.
set -uo pipefail
export PATH=/opt/conda/envs/osmo_refdb/bin:$PATH
# run_annotate_cached.sh <genomes_dir> <strains_dir> <out_dir> <release> [par]: the genomes in
# <genomes_dir>, their cached proteins in <strains_dir>/<bid>/<bid>.faa.gz, the annotation in
# <out_dir>/<bid>/ (in a phenotype_validation directory: genomes/, random_panel/strains/, annotate/).
GENOMES=${1:?usage: run_annotate_cached.sh <genomes_dir> <strains_dir> <out_dir> <release> [par]}
STRAINS=${2:?missing strains_dir}
OUT=${3:?missing out_dir}
REL=${4:?missing release dir}
PAR=${5:-8}
mkdir -p "$OUT"

one() {
  g="$1"; STRAINS="$2"; OUT="$3"; REL="$4"
  b=$(basename "$g"); bid=${b%%_*}
  d="$OUT/$bid"
  [ -s "$d/$bid.gene_counts.tsv" ] && { echo "$bid cached"; return 0; }
  prot="$STRAINS/$bid/$bid.faa.gz"
  [ -s "$prot" ] || { echo "$bid rc=no_cached_proteins"; return 0; }
  mkdir -p "$d"
  t0=$SECONDS
  input="$g"
  case "$g" in
    *.gz) input="$d/$bid.input.fna"
          gunzip -c "$g" > "$input" || { echo "$bid rc=gunzip_failed"; rm -f "$input"; return 0; } ;;
  esac
  gunzip -c "$prot" > "$d/$bid.input.faa" || { echo "$bid rc=protein_gunzip_failed"; rm -f "$d/$bid.input.faa"; return 0; }
  osmotool annotate "$REL" "$input" --proteins "$d/$bid.input.faa" \
    --out_prefix "$d/$bid" --threads 1 >"$d/$bid.log" 2>&1
  rc=$?
  [ "$input" != "$g" ] && rm -f "$input"
  rm -f "$d/$bid.input.faa"
  echo "$bid rc=$rc wall=$((SECONDS-t0))s"
  [ $rc -ne 0 ] && tail -5 "$d/$bid.log"
  return 0
}
export -f one


# CHUNK="i/K" (SLURM arrays): process only the genomes with (position - 1) % K == i,
# so K array tasks split the sorted genome list between them without overlap.
chunk_filter() {
  if [ -n "${CHUNK:-}" ]; then awk -v i="${CHUNK%/*}" -v k="${CHUNK#*/}" '(NR - 1) % k == i'; else cat; fi
}

T0=$SECONDS
INPUTS="$OUT/.inputs${CHUNK:+.${CHUNK%/*}}"
find "$GENOMES/" -maxdepth 1 \( -name '*.fna' -o -name '*.fna.gz' \) | sort | chunk_filter > "$INPUTS"
xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {} "$STRAINS" "$OUT" "$REL" < "$INPUTS"
N=$(wc -l < "$INPUTS" | tr -d ' '); W=$((SECONDS-T0))
echo "=== ${PAR}x1 (cached proteins): ${W}s wall for $N genomes = $(python3 -c "print(f'{$W/$N:.2f}')")s/genome ==="
