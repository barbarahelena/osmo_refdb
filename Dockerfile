FROM continuumio/miniconda3:24.1.2-0

LABEL org.opencontainers.image.description="osmo_refdb: build & benchmark reference DIAMOND/HMM databases for osmoadaptation genes"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

# --- conda environment (mafft, trimal, hmmer, orfm, diamond, cd-hit, wgsim, python deps) ---
COPY environment.yml /opt/osmo_refdb/environment.yml
RUN conda env create -f /opt/osmo_refdb/environment.yml \
    && conda clean -afy

# --- install osmotool itself (needed for the DIAMOND side of the benchmark,
#     via `osmotool profile`). By default this installs from a sibling
#     ../osmotool checkout (for local testing before osmotool is published);
#     once osmotool is on PyPI, replace with `pip install osmotool==<version>`.
COPY osmotool_src /opt/osmotool_src
RUN conda run -n osmo_refdb env SETUPTOOLS_SCM_PRETEND_VERSION=0.1.0 \
    pip install --no-deps /opt/osmotool_src

COPY . /opt/osmo_refdb/

ENV PATH="/opt/conda/envs/osmo_refdb/bin:$PATH"
WORKDIR /opt/osmo_refdb

ENTRYPOINT ["/bin/bash"]
