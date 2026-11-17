#!/bin/bash
# The canonical build of a release: run_pipeline.sh all, from a committed
# checkout, inside the pinned image, followed by PROVENANCE.json. Local runs are
# for developing; this is the run a release is made from.
#
#   source hpc/slurm/env.sh
#   export REUSE_REFS_FROM=releases/v10        # the refs/alignments snapshot to build from
#   sbatch hpc/slurm/05_pipeline.sbatch
#
# From the snapshot, everything after the fetch (split, alignment of changed
# families, HMMs, calibration, DIAMOND, benchmark) is rerun here. A fetch from
# scratch (no REUSE_REFS_FROM) needs the login node's internet and produces a
# different database vintage: UniProt/RefSeq are live; see hpc/README.md.
# Resource requests are placeholders: no figures are known for this cluster yet.
#SBATCH --job-name=osmo_refdb_build
#SBATCH --cpus-per-task=16        # ADJUST
#SBATCH --mem=64G                 # ADJUST
#SBATCH --time=24:00:00           # ADJUST
#SBATCH --output=osmo_refdb_build_%j.log
set -euo pipefail
cd "${REPO}"
: "${RELEASE:?set RELEASE (the new release name)}"
: "${SIF:?set SIF}"
[ -e "releases/${RELEASE}" ] && { echo "ERROR: releases/${RELEASE} exists; use a new name or move it" >&2; exit 1; }

# a build is only attributable to a commit if the tree is clean
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "WARNING: the working tree is dirty; PROVENANCE.json will say so" >&2
fi
echo "commit: $(git rev-parse HEAD)  node: $(hostname)  arch: $(uname -m)  AVX2 cores: $(grep -c avx2 /proc/cpuinfo || true)"
SIF_SHA256=$(sha256sum "${SIF}" | cut -d' ' -f1)

apptainer exec --pwd /work \
    --env "SIF_SHA256=${SIF_SHA256},THREADS=${SLURM_CPUS_PER_TASK},JOBS=${SLURM_CPUS_PER_TASK},REUSE_REFS_FROM=${REUSE_REFS_FROM:-},SPLIT_MODE=${SPLIT_MODE:-random}" \
    "${SIF}" bash -c '
        set -euo pipefail
        bash run_pipeline.sh all "'"${RELEASE}"'"
        python pipeline/write_provenance.py --release "releases/'"${RELEASE}"'"
    '
echo "done: releases/${RELEASE}; compare two builds with: python hpc/verify_reproducible.py releases/A releases/B"
