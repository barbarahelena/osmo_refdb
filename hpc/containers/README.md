# Containers: what each step needs

Each module under `hpc/nextflow/modules/local/` declares its own container, with what the step calls
written above it. They are placeholders (`container 'todo.invalid/<module>'`) until the images exist;
`grep -rn todo.invalid hpc/nextflow/modules` lists what is still open. Versions below are those of the
validated builds (the local arm64 image, and `.venv-nacl` for the analysis), so a new image changes as
little as possible; none of them is a requirement.

| process | calls | image |
|---|---|---|
| `PIN` | git, sha256sum | none: runs on the node |
| `FETCH_GMGC` | python 3 (stdlib only), bash | python 3.11 |
| `FETCH_BACDIVE` | curl, python 3 (stdlib only), awk, gzip, xargs, sed, grep | curl, python 3.11, gawk, gzip, findutils, sed, grep |
| `CATALOG_SCREEN` (`SCREEN`) | python + pyyaml | python 3.11, pyyaml 6.0.3 |
| `PREPARE_BACDIVE` | bash, ls, grep, ln, cp | anything small (the `SCREEN` image will do) |
| `BACDIVE_PROTEINS` | prodigal, hmmsearch, awk, find, xargs, sort, gzip | prodigal 2.6.3, hmmer 3.4, gawk, findutils, gzip. The osmotool image has prodigal and hmmer too, and using it keeps Prodigal identical to the one `osmotool annotate` uses (the script says it must run as `-p meta` exactly like it) |
| `ANALYZE_NACL` | python + numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, xgboost, pyyaml | python 3.11, numpy 2.4.6, pandas 3.1.1, scipy 1.17.1, scikit-learn 1.9.1, matplotlib 3.11.2, xgboost 3.2.0 (>= 3.1: older versions misread the saved model's base score), pyyaml 6.0.3, statsmodels (not pinned: the venv that ran `analyze.py` is not recorded) |
| `INDEX_CATALOG` | `osmotool index-catalog --aligner bwa` | **osmotool image + classic bwa** (see below) |
| `ANNOTATE_CATALOG` | `osmotool build-catalog`: hmmsearch, bwa-mem3; python for `14_package_catalog.py` | **osmotool image** with catalog mode |
| `BACDIVE_ANNOTATE` | `osmotool annotate`: hmmsearch, diamond, prodigal | osmotool image |
| `PROFILE_SAMPLE` | `osmotool profile`: diamond, hmmsearch, orfm, bwa-mem3 | osmotool image |
| `BUILD_RELEASE`, `PROVENANCE` | the whole toolchain **and** osmotool (the benchmarks run `osmotool profile` and `annotate`) | **osmotool image + the build tools** (see below) |

`PROVENANCE` is in the build image because `write_provenance.py` runs every tool's `--version`.

## osmotool is a container, not a package

`barbarahelena/osmotool` on Docker Hub, built by `docker-publish.yml` from every `v*.*.*` tag. It holds
python 3.11, diamond, prodigal, hmmer, orfm, bwa-mem3, pandas and xgboost. Two things are missing for this
pipeline, both in the osmotool repo:

1. **A tag with catalog mode.** The latest tag is v0.8.0. `catalog-profile` (build-catalog, index-catalog,
   `profile --catalog`) is 3 commits ahead of `main`, and the report-unit change (PR #13) is open. Until
   that is merged and tagged, no image has the commands `ANNOTATE_CATALOG`, `INDEX_CATALOG` and
   `PROFILE_SAMPLE` call.
2. **Classic `bwa`.** `environment.yml` has `bwa-mem3` only. `INDEX_CATALOG` builds the full index with
   classic bwa (bwa-mem3 would need 100+ GB for the whole catalog), so `- bwa` has to be added there. The
   repo's own `environment.yml` had the same gap; it is fixed here.

## The build image is the one that Seqera Containers cannot make alone

`BUILD_RELEASE` needs osmotool and about a dozen other tools in one image. Seqera Containers builds from
package names, and osmotool is not a package, so this one needs osmotool's image as the base:

- **A Dockerfile on top of the osmotool image** (simplest; build and push it like osmotool's):
  ```dockerfile
  FROM barbarahelena/osmotool:<tag with catalog mode>
  RUN conda install -y -n osmotool -c conda-forge -c bioconda \
        mafft=7.526 trimal=1.5.1 cd-hit=4.8.1 wgsim=1.0 insilicoseq=2.0.1 \
        scikit-learn=1.9.0 scipy=1.17.1 matplotlib=3.11.0 biopython=1.87 \
        pyyaml=6.0.3 requests=2.34.2 statsmodels gawk \
      && conda clean -afy
  RUN apt-get update && apt-get install -y --no-install-recommends git curl \
      && rm -rf /var/lib/apt/lists/*
  ```
  What osmotool's image already provides (hmmer, diamond, prodigal, orfm, bwa-mem3, bwa once added, pandas,
  xgboost, procps, python 3.11) is not repeated. The validated builds had hmmer 3.4, diamond 2.2.4,
  prodigal 2.6.3, orfm 2.1.1, bwa-mem3 0.17.0, bwa 0.7.19, pandas 3.0.3.
- **Wave** can add conda packages on top of an existing image (its CLI/API, not the Seqera Containers web
  form); check its documentation for the exact call.
- Not recommended now: splitting `BUILD_RELEASE` into one process per pipeline step, so that no process
  needs osmotool and the toolchain together. It would also remove the largest single image, but it means
  taking `run_pipeline.sh all` apart.

## Not verified here

Nothing here was built or run against a registry. bwa-mem3 needs AVX2 on x86-64: the validated image was
arm64, so an x86-64 build is untested. Seqera Containers can build linux/amd64.
