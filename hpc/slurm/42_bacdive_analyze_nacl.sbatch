#!/bin/bash
# The BacDive analysis on the new annotation, then the NaCl model. Needs every task of
# 41 done. Measured on the laptop: build_matrix + analyze + analyze_extended about 16 min
# (v10), NaCl train about 20 min + explain about 10 min.
#SBATCH --job-name=osmo_bacdive_nacl
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G                 # peak memory not recorded; ran on a 7.65 GiB laptop VM
#SBATCH --time=03:00:00
#SBATCH --output=osmo_bacdive_nacl_%j.log
set -euo pipefail
cd "${REPO}"
PV="releases/${RELEASE}/phenotype_validation"
apptainer exec --pwd /work --env "THREADS=${SLURM_CPUS_PER_TASK}" "${SIF}" bash -c "
    set -euo pipefail
    python pipeline/phenotype_validation/build_matrix.py --dir ${PV}
    python pipeline/phenotype_validation/analyze.py --dir ${PV}
    python pipeline/phenotype_validation/analyze_extended.py --dir ${PV}
    if [ -d ${PV}/strategy_control/annotate ]; then
        python pipeline/phenotype_validation/analyze_strategy_control.py --dir ${DATA_DIR} --osmo-dir ${PV}
    fi
    PYTHON=python bash run_pipeline.sh nacl_model ${RELEASE}
    python pipeline/write_provenance.py --release releases/${RELEASE}
"
