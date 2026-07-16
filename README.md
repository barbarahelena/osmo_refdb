# osmo_refdb

Build and benchmark reference DIAMOND + HMM databases for osmoadaptation
gene families (ectA, ectB, ectC, betL, kdpA, nhaA, ...). This is the
**maintainer-facing** tool: use it to (re)build the reference database that
[`osmotool`](https://github.com/barbarahelena/osmotool) consumes for
profiling, and to benchmark DIAMOND vs. HMM-based detection whenever you
add/change a gene family.

This is a separate project from `osmotool` on purpose: `osmotool` is the
lightweight, publication-facing CLI tool that end users run against a
frozen, versioned reference database release. `osmo_refdb` is the
reproducible pipeline that produces that release. You should not need
MAFFT, HMMER, cd-hit, InSilicoSeq, etc. installed just to run `osmotool`.

## Quick start (Docker)

```bash
docker build -t osmo_refdb .
docker run --rm -v "$(pwd):/opt/osmo_refdb" osmo_refdb -c "bash run_pipeline.sh all v1"
```

## Quick start (conda)

```bash
conda env create -f environment.yml
conda activate osmo_refdb
bash run_pipeline.sh all v1
```

## What `run_pipeline.sh` does

```
bash run_pipeline.sh build v1        # fetch, cluster, split, align,
                                      # hmmbuild, calibrate, press, diamond makedb
bash run_pipeline.sh benchmark v1    # simulate reads, run DIAMOND+HMM, score
bash run_pipeline.sh all v1          # build then benchmark
```

Everything for a given release name is written to `releases/<name>/`
(git-ignored), so you can keep multiple versions around and reproduce any
of them independently.

## Adding a new gene family

Edit `families.yaml` — add one entry with a `name`, `positive_query`
(UniProt REST query) and `negative_query` (a hard-negative Pfam-family
query). Then rerun `bash run_pipeline.sh all <release_name>`. No other file
needs to change.

See `families.yaml` for query syntax notes (in particular: UniProt's REST
API does **not** support free-text `family:"..."` queries — use
`xref:pfam-PFxxxxx` instead, verified via the
[InterPro/Pfam API](https://www.ebi.ac.uk/interpro/api/entry/pfam/)).

## Pipeline steps (`pipeline/`)

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_fetch_refs.py` | Fetch positive + hard-negative sequences per family from UniProt |
| 2 | `02_cluster_cdhit.sh` | CD-HIT cluster positives at 90% identity (remove redundancy) |
| 3 | `03_split_train_test.py` | Split positives into train (build) / test (held-out benchmark) |
| 4 | `04_align_trim.sh` | MAFFT align + trimAl trim the TRAIN positives |
| 5 | `05_build_hmms.sh` | hmmbuild per family; score positive-test/negative sets |
| 6 | `06_calibrate_cutoffs.py` | Set per-family HMM GA (gathering) cutoffs |
| 7 | `07_press_hmms.sh` | Concatenate + hmmpress into one binary HMM database |
| 8 | `08_build_diamond_db.sh` | Build a DIAMOND db from TRAIN positives only |
| 9 | `09_simulate_reads.py` | Simulate reads from held-out TEST positives + negatives |
| 10 | `10_run_benchmark.sh` | Run DIAMOND (`osmotool profile`) + HMM (hmmscan) on simulated reads |
| 11 | `11_compute_metrics.py` | Precision/recall/F1, ROC/PR curves, per-family best-threshold summary |

Steps 4–8 only ever touch the TRAIN split; steps 9–11 only ever touch the
TEST split — this train/test separation is required so the benchmark
measures genuine generalization rather than memorization.

## Releasing a database for `osmotool` to use

Once you're happy with a release's benchmark results, publish
`releases/<name>/osmo_refdb.dmnd` and `releases/<name>/hmms/osmo_refdb.hmm`
(+ its `.h3*` index files) as a versioned data release (e.g. a GitHub
release asset or Zenodo record) that `osmotool` users download and point
`--database` at. Cite the specific release version in any publication using
`osmotool`.
