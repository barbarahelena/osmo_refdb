#!/bin/bash
# Release-independent, per genome: size/GC, Prodigal proteins (cached gzipped), and the
# random-panel pool search, for the BacDive genomes and the salt-in control's. A SLURM
# array: task i of NCHUNKS takes every NCHUNKS-th genome. Submit with
#   sbatch --array=0-$((NCHUNKS-1)) hpc/slurm/40_bacdive_proteins.sbatch
# Measured: 19,965 s at 8 parallel jobs for 6,637 genomes = about 44 core-hours in all,
# about 24 core-seconds per genome (about 66 min per task at NCHUNKS=40 on one core).
#SBATCH --job-name=osmo_bacdive_proteins
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G                  # each job ran inside a 7.65 GiB VM with 8 in parallel: well under 1 GB each
#SBATCH --time=02:00:00           # about 17 min expected at 4 cores; margin for slow nodes
#SBATCH --output=osmo_bacdive_proteins_%A_%a.log
set -euo pipefail
cd "${REPO}"
POOL="${DATA_DIR}/random_panel/pool.hmm"
[ -s "${POOL}" ] || { echo "ERROR: ${POOL} missing (snapshot, see hpc/README.md)" >&2; exit 1; }
for PV in "${DATA_DIR}" "${DATA_DIR}/strategy_control"; do
    [ -d "${PV}/genomes" ] || continue
    apptainer exec --pwd /work --env "CHUNK=${SLURM_ARRAY_TASK_ID}/${NCHUNKS}" "${SIF}" \
        bash pipeline/phenotype_validation/run_random_panel.sh "${PV}" "${POOL}" "${SLURM_CPUS_PER_TASK}"
done
