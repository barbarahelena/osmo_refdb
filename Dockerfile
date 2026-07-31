FROM continuumio/miniconda3:24.1.2-0

LABEL org.opencontainers.image.description="osmo_refdb: build & benchmark reference DIAMOND/HMM databases for osmoadaptation genes"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        git \
    && rm -rf /var/lib/apt/lists/*

# --- conda environment (mafft, trimal, hmmer, orfm, diamond, cd-hit, wgsim, python deps) ---
COPY environment.yml /opt/osmo_refdb/environment.yml
RUN conda env create -f /opt/osmo_refdb/environment.yml \
    && conda clean -afy

# --- install osmotool itself (needed for the DIAMOND/HMM side of the
#     benchmark, via `osmotool profile`/`osmotool annotate`). Installed
#     directly from its GitHub repo rather than vendoring a local copy --
#     pinned to a specific commit for reproducibility; bump deliberately
#     when osmo_refdb should pick up newer osmotool changes. Once osmotool
#     has tagged releases, prefer pinning a tag (@v0.1.0) or switching to
#     `pip install osmotool==<version>` once it's on PyPI.
#
# Bumped 4cbc3ba->2826ebc (2026-07-31): the old pin predated ~25 commits of
# real functional changes, including quantifier.KNOWN_FAMILIES being fixed
# from the original ~15-family list to the full v7 43-family panel
# (c097a43) -- mrpG/gshB/trkA/etc. (this branch's own RefSeq/murB test
# subjects) aren't in the old pin's known-family list at all. Confirmed
# this branch's own benchmarks (mrpG/otsA/mscS/otsB, the v9 rebuild, murB's
# PF01565 fix) were run against an ad-hoc environment already close to
# current osmotool HEAD, not the stale pin -- this bump makes the Dockerfile
# match what was actually tested, rather than the other way around.
RUN conda run -n osmo_refdb pip install --no-deps \
    "git+https://github.com/barbarahelena/osmotool.git@2826ebc"

COPY . /opt/osmo_refdb/

ENV PATH="/opt/conda/envs/osmo_refdb/bin:$PATH"
WORKDIR /opt/osmo_refdb

ENTRYPOINT ["/bin/bash"]
