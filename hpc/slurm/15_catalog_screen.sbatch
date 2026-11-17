#!/bin/bash
# The gate before the release is frozen: which families claim catalog genes
# that look like something else (pipeline/15_catalog_screen.py). Read its
# output; a flag means "go and look", and a real confuser class means a decoy
# and a rebuild, before BacDive and the NaCl model are run on these HMMs.
#SBATCH --job-name=osmo_catalog_screen
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G                  # ADJUST
#SBATCH --time=00:30:00
#SBATCH --output=osmo_catalog_screen_%j.log
set -euo pipefail
cd "${REPO}"
apptainer exec --pwd /work "${SIF}" python pipeline/15_catalog_screen.py \
    --release "releases/${RELEASE}" --catalog "${CATALOG_NAME}"
