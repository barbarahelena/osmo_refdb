#!/bin/bash
# Make releases/<release>/phenotype_validation/ point at the release-independent data
# (genomes, cached proteins, random panel), as the BacDive README describes for
# "another release". Run once per release, after the proteins job (40) and before 41.
#   bash hpc/prepare_bacdive_release.sh <release> <data_dir>
set -euo pipefail
REL="${1:?release}"; DATA="$(cd "${2:?data dir}" && pwd)"
PV="releases/${REL}/phenotype_validation"
mkdir -p "${PV}/strategy_control"
for d in genomes random_panel; do ln -sfn "${DATA}/${d}" "${PV}/${d}"; done
for d in genomes random_panel; do [ -e "${DATA}/strategy_control/${d}" ] && ln -sfn "${DATA}/strategy_control/${d}" "${PV}/strategy_control/${d}"; done
for f in pilot_strains.tsv bacdive_classes.tsv; do [ -s "${DATA}/${f}" ] && cp -f "${DATA}/${f}" "${PV}/${f}"; done
[ -s "${DATA}/strategy_control/pilot_strains.tsv" ] && cp -f "${DATA}/strategy_control/pilot_strains.tsv" "${PV}/strategy_control/"
ls -l "${PV}" | tail -n +2
