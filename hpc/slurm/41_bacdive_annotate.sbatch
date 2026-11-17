#!/bin/bash
# osmotool annotate with the release's HMMs, per genome, from the cached Prodigal
# proteins (release-dependent, so this runs after the gate). A SLURM array like 40.
# Measured: 4,365 s at 3 parallel jobs for 6,637 genomes = about 3.6 core-hours in all.
#SBATCH --job-name=osmo_bacdive_annotate
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=osmo_bacdive_annotate_%A_%a.log
set -euo pipefail
cd "${REPO}"
PV="releases/${RELEASE}/phenotype_validation"
for D in "${PV}" "${PV}/strategy_control"; do
    [ -d "${D}/genomes" ] || continue
    apptainer exec --pwd /work --env "CHUNK=${SLURM_ARRAY_TASK_ID}/${NCHUNKS}" "${SIF}" \
        bash pipeline/phenotype_validation/run_annotate_cached.sh "${D}" "releases/${RELEASE}" "${SLURM_CPUS_PER_TASK}"
done
