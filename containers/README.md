# Containers

Every step has its own container, declared in its module the nf-core way. There are no params for it:

```groovy
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/hmmer:3.4--h7d74f8d_5' :
    'biocontainers/hmmer:3.4--h7d74f8d_5' }"
```

Images without a registry in their name come from `quay.io`, which is set in `nextflow.config` as in
nf-core. The osmotool and Seqera images name their registry themselves. There is no `conda`
directive, because osmotool is not on conda. Each module still has an `environment.yml` next to its
`main.nf`, in nf-core format, as the record of what its image holds and at which versions. The
osmotool modules list the osmotool image's own environment, plus osmotool from its git tag.

Each step reports its tool versions to the `versions` topic channel: in its outputs,
`tuple val("${task.process}"), val('<tool>'), eval('<version command>'), topic: versions`. The
workflow collects them from `channel.topic('versions')` into `software_versions.yml`, and
`PROVENANCE` records them per step in `PROVENANCE.json` and copies the file next to the release. This
needs Nextflow >= 25.04. The exception is `OSMOTOOL_INDEXCATALOG`: under storeDir it is skipped once
the index exists, so it writes its versions into the index directory instead.

## Which step runs in which image

| Image | Steps |
|---|---|
| `docker.io/barbarahelena/osmotool:0.8.0` | `OSMOTOOL_ANNOTATE` (aliases `_GENOME`, `_BENCHMARK`), `OSMOTOOL_PROFILE` (and `_RECALL`), `OSMOTOOL_BUILDCATALOG`, `OSMOTOOL_INDEXCATALOG`, `BACDIVE_PROTEINS` and `BACDIVE_ANNOTATE` (Prodigal and hmmsearch must be osmotool's own: the cached proteins are what `annotate --proteins` reads), `BENCHMARK_READS` (osmotool profile, orfm and hmmscan over a whole read set) |
| `biocontainers/orfm:2.1.1--h54198d6_1` | `ORFM` (same build as in the osmotool image) |
| `biocontainers/hmmer:3.4--h7d74f8d_5` | `HMMER_HMMBUILD` (and `_DECOY`), `HMMER_HMMSEARCH`, `HMMER_HMMSCAN` (same build as in the osmotool image) |
| `biocontainers/cd-hit:4.8.1--h5ca1c30_13` | `CDHIT_DECOY` |
| `biocontainers/prodigal:2.6.3--h577a1d6_11` | `POOL_PROTEOMES` (the build in the osmotool image) |
| `biocontainers/wgsim:1.0--h577a1d6_10` | `WGSIM` |
| `biocontainers/python:3.11` | Python standard library only: `CATALOG_FETCH`, `CATALOG_PACKAGE`, `VALIDATION_CHECKSYSTEMS`, `VALIDATION_TARGETCOORDS`, `VALIDATION_READRECALL`, `ALIGNMENTS_RECORD`, `POOL_CANDIDATES`, `POOL_BUILD` |
| Seqera `py-yaml` | `SPLIT_REFS`, `HMM_DECOYSELECT`, `DIAMOND_EXTRAREFS`, `CATALOG_SCREEN`, `BENCHMARK_ANNOTATEPREP`, `BENCHMARK_ANNOTATESCORE` |
| Seqera `refs` | `RELEASE_REFS`, `PFAM_FETCH` |
| Seqera `align` | `MAFFT_ALIGN`, `TRIMAL` (and their `_DECOY` aliases; mafft 7.526 is not on biocontainers) |
| Seqera `hmms` | `HMM_CALIBRATE` (hmmpress, python, numpy, pyyaml; its mafft, trimal and cd-hit are no longer used) |
| Seqera `diamond-py` | `DIAMOND_CALIBRATE` (python, numpy, pyyaml) |
| `biocontainers/diamond:2.2.8--he361c42_0` | `DIAMOND_TRAINREFS`, `DIAMOND_MAKEDB`, `DIAMOND_BLASTP` (the build in the osmotool image) |
| Seqera `simulate` | `BENCHMARK_SIMULATE_CALIB`, `BENCHMARK_SIMULATE_TEST` |
| Seqera `metrics` | `BENCHMARK_METRICS_CALIB`, `BENCHMARK_METRICS_TEST` |
| Seqera `nacl` | `BACDIVE_MATRIX`, `BACDIVE_ANALYZE`, `NACL_MODEL` |
| Seqera `bacdive-fetch` | `BACDIVE_HARVEST`, `BACDIVE_GENOMES`, `POOL_PFAM` |
| Seqera `provenance` | `PROVENANCE` |

## The Seqera images

Made at [seqera.io/containers](https://seqera.io/containers) for linux/amd64, once as Docker (used
by Docker and by Apptainer, which converts it) and once as Singularity (`oras://`, used when the
engine is `singularity`). The packages and versions of each are in the `environment.yml` of the
modules that use it. To change one, edit those files, build the image again there, and replace both
URLs in the module's `container` line.

| Name | Docker | Singularity |
|---|---|---|
| `py-yaml` | `community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81` | `oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580` |
| `refs` | `community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:73160cfec41e18d3` | `oras://community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:ddb8dc66d39a1fc1` |
| `align` | `community.wave.seqera.io/library/trimal_mafft_python_pyyaml:d3ab7e2b4cc1b639` | `oras://community.wave.seqera.io/library/trimal_mafft_python_pyyaml:0aafcb62c9cac1a6` |
| `hmms` | `community.wave.seqera.io/library/cd-hit_hmmer_trimal_mafft_pruned:608e5325648e4136` | `oras://community.wave.seqera.io/library/cd-hit_hmmer_trimal_mafft_pruned:cf68077b56599bf6` |
| `diamond-py` | `community.wave.seqera.io/library/diamond_numpy_python_pyyaml:9cd53c8d486ed198` | `oras://community.wave.seqera.io/library/diamond_numpy_python_pyyaml:5b3b1ed5ce857be9` |
| `simulate` | `community.wave.seqera.io/library/insilicoseq_wgsim_python_pyyaml:e5b214e4c49c8d0d` | `oras://community.wave.seqera.io/library/insilicoseq_wgsim_python_pyyaml:a4bccdbc78a6f5da` |
| `metrics` | `community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scikit-learn:55ad9ba88a0d577d` | `oras://community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scikit-learn:53d0d02acba1bb21` |
| `nacl` | `community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:edc7132220048e3f` | `oras://community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:92ae9dc030c0ee62` |
| `bacdive-fetch` | `community.wave.seqera.io/library/curl_gawk_python:de4685943c4289cd` | `oras://community.wave.seqera.io/library/curl_gawk_python:2d801c4a5165bf61` |
| `provenance` | `community.wave.seqera.io/library/git_python_pyyaml:241a3ffe7fa635b5` | `oras://community.wave.seqera.io/library/git_python_pyyaml:64763cd90b56dd64` |

### Still to do (checked against the registry on 2026-11-17)

- **To rebuild.** Seqera images have no `find`, `xargs`, `gzip`, `zcat` or `cmp`. `bacdive-fetch` needs
  `conda-forge::gzip` and `conda-forge::findutils` (the genome download in `BACDIVE_GENOMES`); its `environment.yml` lists them.
- `python=3.11` resolved to 3.11.0 in every Seqera image. The osmotool image has 3.11.18. To have
  one Python everywhere, pin `python=3.11.18` and rebuild.

## osmotool

`barbarahelena/osmotool:0.8.0` on Docker Hub (digest
`sha256:f6d0cd097390d89446d4975ac9c96fae5b98258f6388d7bfb48b1d058fae8c58`, the same as `latest` on
2026-11-17). Checked here: osmotool 0.8.0, diamond 2.2.8, hmmer 3.4, prodigal 2.6.3, orfm 2.1.1, Python
3.11 with numpy and pandas, and no pyyaml.

- **Catalog mode is not in 0.8.0.** That means `build-catalog`, `index-catalog`, `profile --catalog`,
  and `bwa` / `bwa-mem3`. It is on osmotool's `catalog-profile` branch. `OSMOTOOL_BUILDCATALOG`,
  `OSMOTOOL_INDEXCATALOG` and `-entry profile_samples` (`--catalog`) need a release that has it. Bump
  the tag in the four `modules/local/osmotool/*/main.nf` files.
- **`--no-mount hostfs`.** On a site whose `apptainer.conf` has `mount hostfs = yes` (Snellius
  does), the host's `/opt` covers the image's `/opt/conda`, and every tool is "not found".
  `nextflow.config` passes `--no-mount hostfs`. To check by hand:

```bash
apptainer exec --no-mount hostfs docker://barbarahelena/osmotool:0.8.0 sh -c \
  'osmotool --version; hmmsearch -h | sed -n 2p; diamond version; prodigal -v; orfm --version'
```

## Check

Run once and read the trace (`pipeline_info/*.trace.tsv` has each task's container). `command not
found` means a tool is missing from that step's image. `-stub-run` does not check images.
