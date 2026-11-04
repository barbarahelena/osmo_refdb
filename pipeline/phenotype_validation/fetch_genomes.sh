#!/usr/bin/env bash
# fetch_genomes.sh <outdir> [bacdive_id...]
#
# Download assemblies for strains in <outdir>/pilot_strains.tsv into
# <outdir>/genomes/ as gzipped FASTA (<bacdive_id>_<accession>.fna.gz). With
# no IDs, fetches every row in the table. PAR (default 4) downloads run at
# once -- NCBI is a shared public service, so keep it modest.
#
# BacDive stores INSDC accessions unversioned (GCA_013177455), but NCBI's path
# needs the version and assembly name (GCA_013177455.1_ASM1317745v1). Rather
# than depend on the `datasets` CLI, list the accession's parent FTP directory
# and read the versioned name back out. Where several versions exist, the
# highest is used: a later version is a corrected assembly of the same strain.
#
# Genomes stay gzipped (~29% of the plain size; ~10 GB instead of ~33 GB for
# the full sweep). osmotool annotate can't read gzip, so run_annotate.sh
# unpacks each one to a temp file for the duration of its job. Each download
# is checked with `gzip -t` and written under a temp name first, so an
# interrupted run never leaves a truncated file that looks cached.
#
# Cached files (.fna.gz, or a plain .fna from the pilot) are left alone.
set -uo pipefail
PV=${1:?usage: fetch_genomes.sh <outdir> [bacdive_id...]}; shift || true
PAR=${PAR:-4}
TSV="$PV/pilot_strains.tsv"
[ -s "$TSV" ] || { echo "missing $TSV -- run build_table.py first" >&2; exit 1; }
mkdir -p "$PV/genomes"

# look columns up by name: positions shift whenever build_table.py gains one
# (the pilot's hard-coded $9 silently became nogrowth_at_pct)
pairs=$(awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) col[$i]=i
          if (!col["bacdive_id"] || !col["insdc"]) { print "ERR" > "/dev/stderr"; exit 1 }
          next }
  { print $col["bacdive_id"] "\t" $col["insdc"] }' "$TSV") \
  || { echo "$TSV lacks a bacdive_id or insdc column" >&2; exit 1; }
if [ "$#" -gt 0 ]; then
  want=" $* "
  pairs=$(printf '%s\n' "$pairs" | awk -F'\t' -v w="$want" 'index(w, " " $1 " ")')
fi

one() {
  bid="$1"; acc="$2"
  out="$G/${bid}_${acc}.fna.gz"
  { [ -s "$out" ] || [ -s "$G/${bid}_${acc}.fna" ]; } && { echo "CACHED $bid"; return 0; }
  p=$(echo "${acc#GCA_}" | sed 's/\(...\)\(...\)\(...\)/\1\/\2\/\3/')
  base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/$p/"
  dirs=$(curl -sS -m 60 --retry 3 --retry-delay 2 "$base" \
         | grep -oE "${acc}\.[0-9]+_[^/\"]+" | sort -u -V -r)
  [ -z "$dirs" ] && { echo "FAIL $bid $acc: no FTP dir (suppressed?)"; return 0; }
  tmp="$out.part"
  # newest version first; a suppressed version keeps its folder but only
  # the metadata reports (the FASTA 404s), so fall back to the next one
  for dir in $dirs; do
    if curl -sS -f -m 600 --retry 3 --retry-delay 2 -o "$tmp" \
         "${base}${dir}/${dir}_genomic.fna.gz" 2>/dev/null && gzip -t "$tmp" 2>/dev/null; then
      mv "$tmp" "$out"; echo "OK $bid $dir"; return 0
    fi
    rm -f "$tmp"
  done
  echo "FAIL $bid $acc: no genomic FASTA in any version ($(echo $dirs | tr ' ' ','))"
}
export -f one
export G="$PV/genomes"

# each line is "<bacdive_id>\t<accession>"; -L 1 hands both fields to one()
LOG="$PV/genome_fetch.log"
printf '%s\n' "$pairs" | grep -v '^$' \
  | xargs -P "$PAR" -L 1 bash -c 'one "$1" "$2"' _ | tee -a "$LOG" \
  | awk '{n[$1]++} NR%250==0 {printf "  %d done (OK=%d CACHED=%d FAIL=%d)\n", NR, n["OK"], n["CACHED"], n["FAIL"]; fflush()}
         $1=="FAIL" {print; fflush()}
         END {printf "downloaded=%d cached=%d failed=%d\n", n["OK"], n["CACHED"], n["FAIL"]}'
