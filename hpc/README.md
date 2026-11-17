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

## Steps

1. **Image** (login node with internet; Apptainer >= 1.2 for `--build-arg`):
   ```bash
   git clone -b release-catalog https://github.com/barbarahelena/osmo_refdb.git && cd osmo_refdb
   apptainer build --fakeroot ~/osmo_refdb_hpc.sif hpc/apptainer/osmo_refdb_hpc.def
   ```
   The image holds the conda env (`environment.yml`: mafft, hmmer, diamond,
   bwa-mem3, ...) and osmotool from branch `catalog-profile`. bwa-mem3 needs
   AVX2 on x86-64 (or ARM NEON); the build script prints whether a node has it.
   Compute nodes usually have no internet: everything is installed here.
2. **Release data** (not in git; about 2 GB from `releases/v10/` on the laptop):
   `osmo_refdb.dmnd`, `osmo_refdb.train_refs.faa`, `osmo_refdb.profile_*`,
   `hmms/`, `refs/`, `markers/`, `results/reads/`, `results/diamond/`,
   `osmo_refdb.profile_report_levels.tsv`. Copy to `<clone>/releases/v10/`
   with `rsync -a`.
3. **Catalog** (login or transfer node):
   ```bash
   bash hpc/00_fetch_gmgc.sh /data/gmgc 40000000   # 40 MB smoke-test sample first
   bash hpc/00_fetch_gmgc.sh /data/gmgc            # the whole human-gut file, CRC-checked
   ```
   `fetch_zenodo_member.py` reads one member out of the 19.3 GB zip with HTTP
   range requests, resumably (`<out>.part`).
4. **Smoke test on the sample**, to size the real job: build a catalog from the
   40 MB sample (`.head`, decompress it and cut at a record boundary first, as
   a truncated gzip ends mid-record) and note `sacct` numbers.
5. **Environment file**: copy `hpc/slurm/env.sh.example` to `env.sh`, edit,
   `source` it. Any path you pass to a job (`R1`, `R2`, `OUT`, the catalog)
   must lie under a directory in `APPTAINER_BIND`.
6. **Build the catalog**: `sbatch hpc/slurm/10_build_catalog.sbatch`. It annotates the catalog
   with the release's HMMs and packages the small hits-only index into
   `releases/<release>/catalogs/<name>/`, then builds the full index in
   `$FULL_INDEX_DIR`.
7. **Profile a sample**:
   `sbatch --export=ALL,R1=...,R2=...,OUT=... hpc/slurm/30_profile_sample.sbatch`.

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
