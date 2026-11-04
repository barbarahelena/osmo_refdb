#!/usr/bin/env bash
# fetch_genomes.sh <outdir> [bacdive_id...]
#
# Download assemblies for strains in <outdir>/pilot_strains.tsv into
# <outdir>/genomes/. With no IDs, fetches every row in the table.
#
# BacDive stores INSDC accessions unversioned (GCA_013177455), but NCBI's path
# needs the version and assembly name (GCA_013177455.1_ASM1317745v1). Rather
# than depend on the `datasets` CLI, list the accession's parent FTP directory
# and read the versioned name back out. Cached files are left alone.
set -uo pipefail
PV=${1:?usage: fetch_genomes.sh <outdir> [bacdive_id...]}; shift || true
TSV="$PV/pilot_strains.tsv"
[ -s "$TSV" ] || { echo "missing $TSV -- run build_table.py first" >&2; exit 1; }
mkdir -p "$PV/genomes"

if [ "$#" -gt 0 ]; then ids="$*"; else ids=$(awk -F'\t' 'NR>1 {print $1}' "$TSV"); fi
ok=0; fail=0; cache=0
for bid in $ids; do
  acc=$(awk -F'\t' -v b="$bid" 'NR>1 && $1==b {print $9}' "$TSV")
  [ -z "$acc" ] && { echo "$bid: not in table"; continue; }
  out="$PV/genomes/${bid}_${acc}.fna"
  [ -s "$out" ] && { cache=$((cache+1)); continue; }
  p=$(echo "${acc#GCA_}" | sed 's/\(...\)\(...\)\(...\)/\1\/\2\/\3/')
  base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/$p/"
  dir=$(curl -sS -m 60 --retry 2 "$base" | grep -oE "${acc}\.[0-9]+_[^/\"]+" | sort -u | head -1)
  if [ -z "$dir" ]; then echo "$bid $acc: no FTP dir (suppressed?)"; fail=$((fail+1)); continue; fi
  if curl -sS -m 300 --retry 2 "${base}${dir}/${dir}_genomic.fna.gz" | gunzip > "$out" 2>/dev/null && [ -s "$out" ]; then
    ok=$((ok+1))
    printf "%s %s: %s (%s bp)\n" "$bid" "$acc" "$dir" "$(grep -v '^>' "$out" | tr -d '\n' | wc -c | tr -d ' ')"
  else
    echo "$bid $acc: download failed ($dir)"; rm -f "$out"; fail=$((fail+1))
  fi
done
echo "downloaded=$ok cached=$cache failed=$fail"
