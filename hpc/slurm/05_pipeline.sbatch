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
# Measured: the same structure (v10-betl-decoy: build + benchmarks + catalogue benchmark)
# took about 2.5 h on 10 CPUs inside a 7.65 GiB VM without running out of memory, so
# 16 cores / 32 GB / 8 h is about 1.6x the cores, 4x the memory and 3x the time.
#SBATCH --job-name=osmo_refdb_build
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=osmo_refdb_build_%j.log
set -euo pipefail
cd "${REPO}"
: "${RELEASE:?set RELEASE (the new release name)}"
: "${SIF:?set SIF}"
[ -e "releases/${RELEASE}" ] && { echo "ERROR: releases/${RELEASE} exists; use a new name or move it" >&2; exit 1; }

# a build can only be tied to one commit if every file matches that commit
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "WARNING: some files differ from the last commit (uncommitted changes); PROVENANCE.json will list them" >&2
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
