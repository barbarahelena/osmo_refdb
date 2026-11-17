# Containers

Two images. Every module names its image in a literal `container` directive, the nf-core way
(`container 'barbarahelena/osmotool:<tag>'`); there are no params for it. Today the tag is `todo`
in every module: once the images are pushed, replace it everywhere with the real tag, or better
with `@sha256:<digest>`:

```bash
grep -rln "barbarahelena/osmo" modules | xargs sed -i 's#:todo#:<tag>#'
```

The registry is left to the container engine's configuration (`docker.io` by default), as in nf-core.

| Image | Processes | Why |
|---|---|---|
| `barbarahelena/osmotool` | `FETCH_GMGC`, `BACDIVE_PROTEINS`, `BACDIVE_ANNOTATE`, `ANNOTATE_GENOME`, `ANNOTATE_CATALOG`, `INDEX_CATALOG`, `PROFILE_SAMPLE`, `CHECK_REFERENCE_SYSTEMS`, `TARGET_COORDS`, `REAL_READ_RECALL` | `osmotool` calls, plus scripts that need only the Python standard library (the image has Python) |
| `barbarahelena/osmo-build` | `RELEASE_REFS`, `SPLIT_REFS`, `ALIGN_FAMILIES`, `BUILD_HMMS`, `BUILD_DIAMOND`, `SIMULATE_READS`, `RUN_BENCHMARK`, `COMPUTE_METRICS`, `BENCHMARK_ANNOTATE`, `ASSEMBLE_RELEASE`, `SIMULATE_GENOME_READS`, `ANALYZE_NACL`, `SCREEN`, `FETCH_BACDIVE`, `PROVENANCE` | build and analysis tools, the numpy/pandas/scikit-learn stack, `pyyaml`, `requests`, `curl`, `gawk` |

A plain Python image is not enough for any of these: the Python steps either need the bio tools next
to them or import `yaml`, `requests`, numpy, pandas, scipy, sklearn, statsmodels, matplotlib or xgboost.

## 1. `osmotool`

`barbarahelena/osmotool` on Docker Hub (built from `v*.*.*` tags). The pipeline needs a tag with

- **catalog mode** (`profile --catalog`, `build-catalog`, `index-catalog`); the latest tag may predate it,
- classic **`bwa`** (0.7.19) next to `bwa-mem3` (the full catalog index is too big for bwa-mem3).

Check:

```bash
apptainer exec docker://barbarahelena/osmotool:<tag> sh -c \
  'osmotool profile --help | grep -i catalog; bwa 2>&1 | head -3; hmmsearch -h | head -2; diamond version; prodigal -v; orfm --version; python3 -V'
```

## 2. `osmo-build`

The osmotool image plus the build, benchmark and analysis tools (`containers/osmo-build/Dockerfile`):

```dockerfile
ARG OSMOTOOL=barbarahelena/osmotool:<tag>
FROM ${OSMOTOOL}
USER root
RUN apt-get update && apt-get install -y --no-install-recommends curl gawk git findutils gzip \
 && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir \
      pyyaml==6.0.3 requests==2.34.2 biopython==1.87 pandas==3.0.3 \
      scikit-learn==1.9.0 scipy==1.17.1 statsmodels matplotlib==3.11.0 \
      xgboost insilicoseq==2.0.1
# Not on pip: mafft 7.526, trimal 1.5.1, cd-hit 4.8.1, wgsim 1.0 (bioconda or from source), on PATH.
```

The versions are those the v10 release was validated with. `PROVENANCE` records `--version` of every
tool, so all of these must be on `PATH`: `osmotool hmmbuild hmmsearch hmmpress diamond prodigal orfm
bwa bwa-mem3 mafft trimal cd-hit wgsim iss git python3`.

## Build and check

```bash
docker build -t barbarahelena/osmo-build:<tag> containers/osmo-build && docker push barbarahelena/osmo-build:<tag>
```

Run once with `-with-trace`: `command not found` in a task means a tool is missing from its image.
`-stub-run` does not check that.

## Not verified

Neither image has been built or run here. `bwa-mem3` needs AVX2 on x86-64; the earlier validated image
was arm64. Pin `xgboost` and `statsmodels` to what `ANALYZE_NACL` was validated with.
