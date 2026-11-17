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
6. **Submit the chain**:
   ```bash
   export RELEASE=v10-canonical REUSE_REFS_FROM=releases/v10
   bash hpc/slurm/submit_all.sh
   ```
   It submits the full index, the release build, the catalog annotation (after
   the build) and the screen (after that), with SLURM dependencies.
7. **The gate**: read the screen's log and `catalogs/<name>/screen.tsv`.
8. **After the gate, on the frozen release**: BacDive annotation as a SLURM
   array and the NaCl retrain (scripts not written yet), the catalog benchmark
   with the real catalog as background (`11d` not extended yet), and
   `sbatch hpc/slurm/30_profile_sample.sbatch` per sample.

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
