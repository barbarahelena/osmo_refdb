# Running the catalog pipeline on an HPC cluster (SLURM + Apptainer)

For building and benchmarking `osmotool profile --catalog` against a full-size
gut gene catalog. Everything here is **untested on a cluster**: the scripts
are small, parametrised and meant to be adjusted on first contact (partition,
memory, time, bind paths).

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
  `pipeline/15_catalog_screen.py`, before BacDive or the NaCl model are run.
- **The gate** is a person reading the screen. A flagged family means "look at
  those genes" (closest relative among the family's own references and the
  suspected confusers). A real confuser class means a decoy and a rebuild; a
  false alarm (the flag also fires on divergent true members) means carry on.
  Only then is the release frozen, and BacDive and the NaCl model run on it.

## Steps

1. **Image** (login node with internet; Apptainer >= 1.2 for `--build-arg`):
   ```bash
   git clone -b <branch> https://github.com/barbarahelena/osmo_refdb.git && cd osmo_refdb
   apptainer build --fakeroot ~/osmo_refdb_hpc.sif hpc/apptainer/osmo_refdb_hpc.def
   ```
   The image holds the conda env (`environment.yml`: mafft, hmmer, diamond,
   bwa-mem3, ...) and osmotool from `OSMOTOOL_REF` (a tag or commit for a real
   release). bwa-mem3 needs AVX2 on x86-64 (or ARM NEON). Compute nodes
   usually have no internet: everything is installed here.
2. **Snapshot** (not in git; about 2 GB from `releases/v10/` on the laptop):
   `refs/`, `alignments/`, and what the build reads from the source release
   (`osmo_refdb.train_refs.faa`, `markers/`, ...). `rsync -a` it to
   `<clone>/releases/<snapshot>/`.
3. **Catalog** (login or transfer node):
   ```bash
   bash hpc/00_fetch_gmgc.sh /data/gmgc 40000000   # 40 MB smoke-test sample first
   bash hpc/00_fetch_gmgc.sh /data/gmgc            # the whole human-gut file, CRC-checked
   ```
   `fetch_zenodo_member.py` reads one member out of the 19.3 GB zip with HTTP
   range requests, resumably (`<out>.part`).
4. **Smoke test on the sample**, to size the real jobs: decompress the `.head`
   file and cut it at a record boundary (a truncated gzip ends mid-record),
   build a catalog from it and note the `sacct` numbers.
5. **Environment file**: copy `hpc/slurm/env.sh.example` to `env.sh`, edit,
   `source` it. Any path you pass to a job (`R1`, `R2`, `OUT`, the catalog)
   must lie under a directory in `APPTAINER_BIND`.
6. **Launch** (everything below is one script, `hpc/launch.sh`):
   ```bash
   source hpc/slurm/env.sh
   export RELEASE=v10-canonical REUSE_REFS_FROM=releases/v10
   bash hpc/launch.sh fetch      # login node: GMGC + BacDive strains and genomes
   bash hpc/launch.sh build      # submits the whole chain up to the gate
   bash hpc/launch.sh analyze    # after the gate: BacDive annotation -> analysis + NaCl model
   ```
   `build` submits, with SLURM dependencies: the full catalog index and the BacDive
   proteins (independent of the release, start at once) and the release build ->
   catalog annotation -> screen.
7. **The gate**: read the screen's log and `catalogs/<name>/screen.tsv`.
8. **Still to do on the frozen release** (not scripted yet): the catalog benchmark with
   the real catalog as background (`11d` is not extended), and per-sample
   `sbatch hpc/slurm/30_profile_sample.sbatch`.

## What is downloaded by a script, and what is not

| input | scripted? | how |
|---|---|---|
| GMGC human gut catalog (2.99 GB) | yes | `launch.sh fetch` -> `00_fetch_gmgc.sh`, resumable, CRC-checked |
| BacDive strain records | yes, no credentials | `00_fetch_bacdive.sh` -> `harvest_bacdive.sh sweep`; batches already present are skipped, so a shipped `bacdive_raw/` snapshot reproduces the published table (BacDive itself is live) |
| BacDive genomes (about 9.4 GB), salt-in control genomes (0.7 GB) | yes | `fetch_genomes.sh` from NCBI FTP, resumable |
| Pfam-A 39.0 | yes, added | `00_fetch_pfam.sh` (pinned URL); only needed to re-derive the random-panel pool |
| random-panel pool (`pool.hmm`, 13 MB) | **no: a snapshot** | re-deriving it needs `screen_proteomes.faa`, which no script creates |
| release snapshot (`refs/`, `alignments/`, ...) | **no: rsync** | UniProt/RefSeq are live, so the snapshot is the input |
| reference genomes for real-genome validation | in git | `test_genomes/` (3), osmotool's `test_genomes/` (9) |
| UniProt / RefSeq (build step 1) | yes, but not used | only for a from-scratch fetch, which is a new database vintage |

## Resources

Every figure is either **measured** (on the laptop: 10 CPUs, 7.65 GiB VM) or **scaled** from a
measured constant by the size of the full catalog (75-86x the 133 Mbp sample, from two
independent gene counts: 13.8 M genes in the Cayman paper, 11.9 M from the file's compressed
size); the basis is stated. Nothing here is a guess about the cluster. Requests carry a
margin; run once and tighten them from `sacct`.

| job | cores | memory | wall | basis |
|---|---|---|---|---|
| fetch GMGC (login node) | 1 | - | about 20-35 min | download rate seen here 1.4-2.7 MB/s for 2.99 GB |
| fetch BacDive + genomes (login node) | 8 | - | about 1.5 h + 50 min | **measured** (harvest 1.5 h; 9.4 GB genome fetch at 8 parallel) |
| release build + benchmarks (`05`) | 16 | 32 GB | request 8 h | **measured** 2.5 h on 10 CPUs inside 7.65 GiB |
| catalog annotation (`10`) | 16 | 24 GB | request 6 h | **scaled**: translation 10-12 min (1 thread), hmmsearch 21-25 min at 4 threads, Python 2.4-2.8 GB |
| full catalog index (`11`, classic bwa) | 1 | 32 GB | request 6 h | **scaled**: classic bwa 92 s and 203 MB for 133 Mbp -> about 2 h and 15-17 GB |
| catalog screen (`15`) | 1 | 8 GB | 30 min | reads one TSV |
| BacDive proteins (`40`) | 40 x 4 | 8 GB each | 2 h each | **measured** 44 core-hours in all (19,965 s at 8 jobs) |
| BacDive annotation (`41`) | 40 x 2 | 4 GB each | 1 h each | **measured** 3.6 core-hours in all (4,365 s at 3 jobs) |
| analysis + NaCl model (`42`) | 8 | 16 GB | request 3 h | **measured** about 16 min analysis + about 30 min NaCl |
| profile, one sample of <= 20M reads (`30`) | 8 | 32 GB | request 6 h | memory **scaled** (mapping loads the index, about 18-21 GB); time **not measured** at this scale |

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

For every job: `sacct -j <id> --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem`.
At minimum:

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
  image can change tool versions. The image bakes in its exact package list
  (`/opt/osmo_refdb/conda-explicit.txt`, `pip-freeze.txt`), is built once per
  release and kept, and its SHA-256 goes into the release's provenance.
- *Code.* Build from a committed checkout, with osmotool installed at a tag or
  commit in the image. Local builds used osmotool from a working tree on
  `PYTHONPATH`, which shadows the installed package; `write_provenance.py`
  detects and records that.
- *Record.* `PROVENANCE.json`, written next to the release by
  `pipeline/write_provenance.py`, is the build's receipt. It changes nothing
  and nothing reads it to build anything; it answers "what exactly made this
  release?" so that two builds can be compared and a difference traced to a
  cause. It holds, in four blocks: **code** (the osmo_refdb commit and any
  uncommitted files, the osmotool commit and whether a working tree shadowed
  the installed one), **environment** (machine type, container and the image's
  checksum, the version of every tool called), **inputs** (families.yaml's
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
  scratch makes a new release, and its `refs/` is then archived. Compute
  nodes may not have internet, so the fetch runs on a login node anyway.
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
1. Build the image once (`hpc/apptainer/osmo_refdb_hpc.def`, `OSMOTOOL_REF` = tag
   or commit) and keep it.
2. Clean checkout of the commit; put the snapshot at `releases/<snapshot>/`.
3. `sbatch hpc/slurm/05_pipeline.sbatch` with `RELEASE` and `REUSE_REFS_FROM`
   set (see the script).
4. Before releasing, rebuild once more from the first build's `refs/` and
   `alignments/` and run
   `python hpc/verify_reproducible.py releases/A releases/B`. It must print
   `REPRODUCED`: HMMs (without DATE/COM lines), cutoffs, cascade config,
   report levels and the benchmark tables identical. The alignments of
   changed families are the one place two *independent* builds can differ,
   which is why the second build reuses the first one's.

## Known gaps

- The catalog benchmark (`pipeline/11d_benchmark_catalog.py`) builds its catalogs
  from the release's own sequences. Adding a real catalog as background (scenarios
  A and B) is the next change; the test reads are back-translated proteins and will
  not map to real genes, so the real catalog tests background competition and
  scale, not real-read mapping.
- `build-catalog` translates every gene in frame 1. 2.1% of the sampled genes have
  a length that is not a multiple of 3 (partial genes) and may be mistranslated; not
  fixed yet.
- `build-catalog` translates in pure Python, single-threaded; the HMM search is
  multi-threaded. Time per stage is not measured at full size.
