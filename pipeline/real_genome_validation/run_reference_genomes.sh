#!/usr/bin/env bash
# run_reference_genomes.sh <release_dir> <genomes_dir> [n_parallel]
#
# Run `osmotool annotate` (default HMM method, Prodigal -p meta) on every
# <genomes_dir>/*.fna and write to
# <release_dir>/real_genome_validation/reference_genomes/<genome>/.
# The nine reference genomes live in osmotool's test_genomes/ (accessions in
# its README). Then compare with expected_systems.tsv:
#   python check_reference_systems.py --release <release_dir>
#
# Uses whichever osmotool is on PATH, or the source tree in OSMOTOOL_SRC
# (python -m osmotool) if that is set. The version is recorded in the output.
set -uo pipefail
REL=${1:?usage: run_reference_genomes.sh <release_dir> <genomes_dir> [par]}
GEN=${2:?missing genomes dir}
PAR=${3:-3}
OUT="$REL/real_genome_validation/reference_genomes"
mkdir -p "$OUT"
if [ -n "${OSMOTOOL_SRC:-}" ]; then
  export PYTHONPATH="$OSMOTOOL_SRC${PYTHONPATH:+:$PYTHONPATH}"
  OSMO="python -m osmotool"
else
  OSMO="osmotool"
fi
export OSMO REL OUT
$OSMO --version > "$OUT/osmotool_version.txt" 2>&1
[ -n "${OSMOTOOL_SRC:-}" ] && git -C "$OSMOTOOL_SRC" describe --tags --always >> "$OUT/osmotool_version.txt" 2>/dev/null
one() {
  g="$1"; name=$(basename "$g" .fna); d="$OUT/$name"; mkdir -p "$d"
  $OSMO annotate "$REL" "$g" --out_prefix "$d/$name" --threads 2 > "$d/$name.log" 2>&1
  echo "$name rc=$?"
}
export -f one
ls "$GEN"/*.fna | xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {}
