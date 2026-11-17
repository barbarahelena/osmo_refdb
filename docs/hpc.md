# Running the whole pipeline on an HPC cluster (Nextflow, SLURM + Apptainer)

For the canonical release build, the catalog, the BacDive validation and the NaCl model, and for
benchmarking `osmotool profile --catalog` against a full-size gut gene catalog. The workflow is at
the repo root (nf-core layout: `main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/`, `conf/`).
It has been run in `-stub-run` mode only (last: Nextflow 25.10.3, after the move to one container per
step), so the wiring (what runs when, what waits for what, the gate, storeDir) is checked; **nothing has
run on a cluster**, and the Seqera images of the build and analysis steps are still to be made
(`containers/README.md`). Resources, bind paths and queues are meant to
be adjusted on first contact.

## What is measured, and what is not

Measured on this repo's laptop setup, on a **40 MB sample of the real human-gut
catalog** (the first 159,560 genes of `GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz`,
137 MB of FASTA, v10 HMMs, 2 threads, while another job was running):

| | measured |
|---|---|
| `build-catalog` (translate, hmmsearch, bwa-mem3 index), wall time | 108 s |
| genes with any HMM hit | 7,889 of 159,560 (4.9%) |
| of which assigned to a panel family | 691 (0.43%) |
| of which sub-cutoff homologs / decoy-model genes | 5,194 / 2,004 |
| full bwa-mem3 index | 472 MB (3.4 bytes per FASTA byte) |
| format | no duplicate IDs, mean gene length 832 nt, 2.1% of genes not a multiple of 3 |

- The sample is the *start* of the file; if the file is ordered in any way, it
  may not be representative of the whole catalog.
- **Not measured**: index build memory and time for the whole catalog, mapping
  speed and peak memory against it, annotation time for the whole catalog.
  These are what the first cluster runs are for. Nothing here extrapolates
  them.
- Real numbers about the source: the human-gut file is 2.99 GB (gzip); the
  Zenodo zip it sits in is 19.3 GB.

## Order of work, and what depends on what

The HMMs decide everything downstream, so they are the first thing built and
the last thing that may still change:

```
release build (HMMs, DIAMOND db, benchmarks)
   |-- catalog annotation -- screen -- GATE -- freeze
   |-- BacDive annotation -- NaCl model
   `-- catalog / profile / annotate benchmarks
full catalog index   (needs only the gene file: starts at once, in parallel)
```

- The catalog and BacDive do not feed each other, but both are made *with the
  HMMs*: a change to the HMMs (say a new decoy) repeats both. So the one that
  can change the HMMs comes first. The catalog annotation does: it is the
  widest real-data test of the HMMs' specificity there is. On 159,560 real gut
  genes it showed gor claiming 127 genes (18% of everything annotated), 86% of
  them scoring below the weakest true gor protein: lipoamide dehydrogenase and
  relatives, which the benchmark never contained. That is found in minutes by
  `modules/local/catalog/screen/templates/15_catalog_screen.py`, before BacDive or the NaCl model are run.
- **The gate** is a person reading the screen. A flagged family means "look at
  those genes" (closest relative among the family's own references and the
  suspected confusers). A real confuser class means a decoy and a rebuild; a
  false alarm (the flag also fires on divergent true members) means carry on.
  Only then is the release frozen, and BacDive and the NaCl model run on it.

## Steps

1. **Containers.** One per step, declared in its module (nf-core style): the osmotool image, biocontainers
   and Seqera images; `containers/README.md` lists which step uses which. The Seqera ones are
   `TODO-seqera-...` placeholders until they are made. Apptainer runs with `--no-mount hostfs`
   (`nextflow.config`): a site that mounts every host filesystem would hide the osmotool image's
   `/opt/conda`. So every directory outside the repo a step reads or writes goes in `extra_binds`.
2. **Snapshot** (not in git; about 2 GB from `releases/v10/` on the laptop): `refs/`, `alignments/`, and
   what the build reads from the source release (`osmo_refdb.train_refs.faa`, `markers/`, ...). `rsync -a`
   it to `<clone>/releases/<snapshot>/`. Also `<repo>/data/phenotype_validation/random_panel/pool.hmm`
   (+ `pool_families.tsv`, `Pfam.version`) and, to reproduce the published BacDive table rather than
   today's, `bacdive_raw/`.
3. **Settings**: edit the `params` block in `nextflow.config` (`catalog_genes`, `full_index_dir`, the
   images, and `queue`/`slurm_extra`/`work_dir`). Paths outside the clone that jobs use go in
   `extra_binds`. The release name is given per run with `--release`.
4. **Run** (from the clone's root, inside a small SLURM job of your own that only runs Nextflow, so nothing
   runs on a login node):
   ```bash
   nextflow run main.nf --release v11
   # at a site with an nf-core/configs profile:  ... --release v11 -profile <institution>
   ```
   It downloads the GMGC catalog and BacDive, then runs in parallel: the full catalog index, the BacDive
   proteins (an array of `nchunks`), and the release build -> catalog annotation -> screen. It stops there.
5. **The gate**: the screen's report is printed at the end of the head job's log, with
   `releases/<release>/catalogs/<name>/screen.tsv` to read. A flag means "look at those genes"; a real
   confuser class means a decoy and a rebuild (a new `release` name); a false alarm means carry on:
   ```bash
   nextflow run main.nf --release v11 -resume --gate_passed true
   ```
   The cached steps are skipped; BacDive annotation (an array), the analysis, the NaCl model and the
   provenance receipt follow.
6. **Samples**: `-entry profile_samples --samples samples.csv` (columns `id,r1,r2`) profiles each sample
   as its own job against the release's catalog and the full index; results in `<outdir>/profile`.
7. **Still to do on the frozen release** (not in the workflow): the catalog benchmark with the real
   catalog as background (`11d` is not extended).

**Skipping finished work: `storeDir` and `-resume`.** Nextflow has two ways to avoid redoing work.

- `-resume` reuses a task's result from `work/` if its inputs and code are unchanged. It is lost if
  `work/` is deleted.
- `storeDir` stores a task's output in a permanent directory of your choosing. If the output is already
  there, the task is not run at all, with or without `-resume` and even if `work/` is gone.

`storeDir` is used where a task's output is a real file that is expensive, does not depend on the
release, and is written once: `CATALOG_FETCH` (the gene file, next to `catalog_genes`) and `OSMOTOOL_INDEXCATALOG`
(the full catalog index, in `full_index_dir`).

Everything else goes from step to step as outputs and inputs; no step reads a published file or
writes into the clone. Results are published (`publishDir`, `conf/modules.config`):

- the release into `releases/<release>/`, by the steps that make its parts: refs, alignments and HMMs,
  the DIAMOND files, the cascade config and report levels, the benchmark results, `catalogs/<name>/`,
  `real_genome_validation/`, `phenotype_validation/`, `nacl_model/`, `PROVENANCE.json`;
- the BacDive tables and harvested batches into `<repo>/<data_dir>` (not the genomes, which stay in
  `work/` for `-resume`).

The one directory osmotool reads (its `DATABASE`) is `RELEASE_DB`: the DIAMOND files, `hmms/` and the
cascade, passed to every osmotool step. As with any output directory, publishing into a release name
that exists overwrites its files: name a new release for every build. The BacDive genomes and proteins are only reused through
`-resume`: deleting `work/` means downloading them again, unless `--bacdive_snapshot` points at
harvested batches (then only the genomes are fetched again).

`pipeline_info/` gets a trace per invocation, including each task's image, peak memory and time: use
it to tighten `conf/base.config`.

## What is downloaded by a script, and what is not

| input | scripted? | how |
|---|---|---|
| GMGC human gut catalog (2.99 GB) | yes | `CATALOG_FETCH` (storeDir) -> `modules/local/catalog/fetch/templates/fetch_zenodo_member.py`, resumable, CRC-checked |
| BacDive strain records | yes, no credentials | `BACDIVE_FETCH` (its template) -> `harvest_bacdive.sh sweep`; with `--bacdive_snapshot <dir>/bacdive_raw` its batches are not fetched again, so a shipped snapshot reproduces the published table (BacDive itself is live). The batches are published to `<repo>/<data_dir>/bacdive_raw` |
| BacDive genomes (about 9.4 GB), salt-in control genomes (0.7 GB) | yes | `fetch_genomes.sh` from NCBI FTP, resumable |
| Pfam-A 39.0 | yes, added | `00_fetch_pfam.sh` (pinned URL); only needed to re-derive the random-panel pool |
| random-panel pool (`pool.hmm`, `pool_families.tsv`, 13 MB) | **no: a snapshot** (`--random_panel`, default `<repo>/<data_dir>/random_panel`) | re-deriving it needs `screen_proteomes.faa`, which no script creates |
| release snapshot (`refs/`, `alignments/`, ...) | **no: rsync** | UniProt/RefSeq are live, so the snapshot is the input |
| reference genomes for real-genome validation | in git | `test_genomes/` (3), osmotool's `test_genomes/` (9); validated by the workflow for every release |
| UniProt / RefSeq (build step 1) | yes, but not used | only for a from-scratch fetch, which is a new database vintage |

## Resources

Every figure is either **measured** (on the laptop: 10 CPUs, 7.65 GiB VM) or **scaled** from a
measured constant by the size of the full catalog (75-86x the 133 Mbp sample, from two
independent gene counts: 13.8 M genes in the Cayman paper, 11.9 M from the file's compressed
size); the basis is stated. Nothing here is a guess about the cluster. Requests carry a
margin; run once and tighten them from `pipeline_info/*.trace.tsv`.

| job | cores | memory | wall | basis |
|---|---|---|---|---|
| `CATALOG_FETCH` | 1 | - | about 20-35 min | download rate seen here 1.4-2.7 MB/s for 2.99 GB |
| `BACDIVE_FETCH` | 8 | - | about 1.5 h + 50 min | **measured** (harvest 1.5 h; 9.4 GB genome fetch at 8 parallel) |
| `BUILD_RELEASE` | 16 | 32 GB | request 8 h | **measured** 2.5 h on 10 CPUs inside 7.65 GiB |
| `OSMOTOOL_BUILDCATALOG` + `CATALOG_PACKAGE` | 16 | 24 GB | request 6 h | **scaled**: translation 10-12 min (1 thread), hmmsearch 21-25 min at 4 threads, Python 2.4-2.8 GB |
| `OSMOTOOL_INDEXCATALOG` (classic bwa) | 1 | 32 GB | request 6 h | **scaled**: classic bwa 92 s and 203 MB for 133 Mbp -> about 2 h and 15-17 GB |
| `CATALOG_SCREEN` | 1 | 8 GB | 30 min | reads one TSV |
| `BACDIVE_PROTEINS` | 40 x 4 | 8 GB each | 2 h each | **measured** 44 core-hours in all (19,965 s at 8 jobs) |
| `BACDIVE_ANNOTATE` | 40 x 2 | 4 GB each | 1 h each | **measured** 3.6 core-hours in all (4,365 s at 3 jobs) |
| `BACDIVE_ANALYZE` + `NACL_MODEL` | 8 | 16 GB | request 3 h | **measured** about 16 min analysis + about 30 min NaCl |
| `OSMOTOOL_PROFILE`, one sample of <= 20M reads | 8 | 32 GB | request 6 h | memory **scaled** (mapping loads the index, about 18-21 GB); time **not measured** at this scale |

**Why classic bwa for the full index.** On the same 133 Mbp: bwa-mem3 builds in 33 s but peaks at
1.37 GB (10.3 bytes per base), classic bwa takes 92 s and peaks at 0.2 GB (1.5 bytes per base).
Scaled linearly to 10-11.5 Gbp that is about 100-120 GB for bwa-mem3 (and its own documentation
gives a higher per-base cost above about 1 Gbp: 72 GB for the human genome, so well over 200 GB)
against about 15-17 GB for classic bwa. Cayman builds its index with classic bwa too. The small
hits-only index shipped in the release stays bwa-mem3.

**Disk** (about 65 GB; ask for 100 GB): GMGC 3 GB, BacDive data 17.4 GB (**measured**), full catalog index
about 18-21 GB (scaled), the release snapshot 2 GB, one build about 8 GB (**measured**: v10 is 7.6 GB),
temporary protein FASTA for the annotation about 4 GB (scaled), the hits-only index and annotation a few GB.

**Total compute**, summed over the table: roughly 44 + 4 core-hours for the BacDive arrays, about 40 for
the release build, about 6 for the catalog annotation, 2 for the index, 5 for the analysis: on the order
of 100 core-hours, none of it needing more than 32 GB per job.

## What to record

Every invocation writes `pipeline_info/<time>_<job>.trace.tsv` (wall time, peak RSS, image, per task),
plus a timeline and a report. At minimum look at:

| job | threads | wall | peak RSS | notes |
|---|---|---|---|---|
| build: hmmsearch + hits index | | | | |
| `index-catalog` (full index) | | | | index size on disk: |
| profile, 20M-read sample, catalog + full index | | | | |
| profile, same sample, no catalog | | | | |
| profile, same sample, catalog + hits index only | | | | |

## Reproducibility: what a canonical HPC build fixes, and what it does not

A release should come from one recorded run, not from whatever a laptop had
that day. Running on the cluster does not make a build reproducible by
itself; these are the things that do, and the things that stay out of reach.

**Fixed by building canonically on the cluster**
- *Architecture.* The local builds so far ran in an arm64 Docker VM
  (`aarch64`); users run x86-64. HMMER and bwa-mem3 use different SIMD code on
  each. Whether scores differ in the last digits has not been measured; a
  canonical x86-64 build removes the question for the release.
- *Environment.* `environment.yml` pins only lower bounds, so rebuilding an
  image can change tool versions. Each module names its image in the repo, so the commit recorded
  in `PROVENANCE.json` fixes which image every step used (and `pipeline_info/` records them per task).
  Use images pinned by a content-addressed tag or digest, and do not rebuild them for a rerun.
- *Code.* Build from a committed checkout (the launch records the commit and a checksum of any
  uncommitted change, and a change to either makes `-resume` rebuild), with osmotool from the image's
  tag. Local builds used osmotool from a working tree on `PYTHONPATH`, which shadows the installed
  package; `write_provenance.py` detects and records that.
- *Record.* `PROVENANCE.json`, written next to the release by
  `modules/local/provenance/templates/write_provenance.py`, is the build's receipt. It changes nothing
  and nothing reads it to build anything; it answers "what exactly made this
  release?" so that two builds can be compared and a difference traced to a
  cause. It holds, in four blocks: **code** (the osmo_refdb commit and any
  uncommitted files, the osmotool commit and whether a working tree shadowed
  the installed one), **environment** (machine type, the version of every tool
  called, in the build image), **inputs** (families.yaml's
  checksum, the UniProt release and fetch dates, which snapshot it was built
  from, a checksum over `refs/` and `alignments/`) and **outputs** (a checksum
  of every HMM, the DIAMOND database, the cutoff files and each catalog's
  annotation and index). Whether two builds agree is checked separately, on the
  files themselves, by `verify_reproducible.py`.

**Not fixed by where it runs**
- *UniProt and RefSeq are live.* Re-running step 1 on the cluster fetches
  whatever they hold that day: a different database vintage, not v10 again.
  The input of a build is therefore the release's `refs/` snapshot
  (`REUSE_REFS_FROM`); the snapshot is shipped with the release. A fetch from
  scratch makes a new release, and its `refs/` is then archived.
- *MAFFT with several threads is not reproducible run to run* (18 of 44 v9
  models differed when realigned from identical input). Alignments are
  therefore inputs: reused from the snapshot, realigned only for families
  whose training set changed (`run_pipeline.sh` does this per family), and
  archived with the release. `hmmbuild` itself is deterministic.

**Measured so far.** Same machine, same refs and reused alignments (v10 vs
v10-betl-decoy): 45 of 46 HMMs identical apart from the `DATE` and `COM` lines
HMMER writes into every file, and the DIAMOND database byte-identical. The one
HMM that differs is betL, which was changed on purpose. **Not measured:** the
same build on arm64 versus x86-64, and the same build twice on the cluster.

**Procedure**
1. Fix the images once (`containers/README.md`: every module's `container` line), commit, and keep them.
2. Clean checkout of the commit; put the snapshot at `releases/<snapshot>/`.
3. `nextflow run main.nf --release <name>` with `reuse_refs_from`
   set; from an empty `work/` (use `-resume` only for the run after the gate).
4. Before releasing, rebuild once more from the first build's `refs/` and
   `alignments/` and run
   `python modules/local/provenance/templates/verify_reproducible.py releases/A releases/B`. It must print
   `REPRODUCED`: HMMs (without DATE/COM lines), cutoffs, cascade config,
   report levels and the benchmark tables identical. The alignments of
   changed families are the one place two *independent* builds can differ,
   which is why the second build reuses the first one's.

## Known gaps

- The catalog benchmark (`modules/local/catalog/package/templates/11d_benchmark_catalog.py`) builds its catalogs
  from the release's own sequences. Adding a real catalog as background (scenarios
  A and B) is the next change; the test reads are back-translated proteins and will
  not map to real genes, so the real catalog tests background competition and
  scale, not real-read mapping.
- `build-catalog` translates every gene in frame 1. 2.1% of the sampled genes have
  a length that is not a multiple of 3 (partial genes) and may be mistranslated; not
  fixed yet.
- `build-catalog` translates in pure Python, single-threaded; the HMM search is
  multi-threaded. Time per stage is not measured at full size.
