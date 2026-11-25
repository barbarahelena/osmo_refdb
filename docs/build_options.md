# Build options

Every option is passed as `--<option>` to `nextflow run main.nf`. The defaults are in `nextflow.config` and the descriptions in `nextflow_schema.json`. Each option that changes the release is recorded in its `PROVENANCE.json`.

```bash
nextflow run main.nf --outdir /data/osmo_refdb/v11 --release v11 --refs_store /project/osmo/refs_store
```

## Input and output

| Option | Default | What it does |
|---|---|---|
| `outdir` | required | The release folder (absolute path, on shared storage) |
| `release` | required | The release's name, written into `PROVENANCE.json`, the catalog and the NaCl model |
| `families_file` | `families.yaml` | The gene families to build; point it at a subset to test-drive new families |
| `markers_file` | `markers.yaml` | The marker genes for `markers/` |
| `catalog_genes` | required | The gene catalog FASTA (absolute path); downloaded there if missing |
| `data_dir` | `bacdive_data` | Folder under `outdir` for the release-independent BacDive data |
| `publish_dir_mode` | `copy` | How results are published (Nextflow `publishDir` mode) |
| `work_dir` | `./work` | Nextflow's work directory; put it on scratch |

## References and split

| Option | Default | What it does |
|---|---|---|
| `split_mode` | `random` | `random`: shuffle sequences into train, calib and test. `taxonomy`: hold out whole genera, to test detection of divergent homologs |
| `refs_store` | none | Keep the UniProt fetch here and reuse it in later builds within the same UniProt release |
| `random_background_n` | 5000 | Proteins in the benchmark's random bacterial background |

## Benchmark

| Option | Default | What it does |
|---|---|---|
| `read_lengths` | none | Also simulate the test reads at these exact lengths, e.g. `100,150,250,300`; scored per length in `results/metrics/summary_by_read_length.tsv`. The calibration reads keep the default length, so the cascade and report levels do not change |
| `background_r1`, `background_r2` | none | A real metagenome with no known osmoadaptation genes, added to the test reads as background; scored in `results/metrics/background_fpr.tsv` |
| `diamond_min_subject_cover` | 0.0 (off) | osmotool profile's `--min_subject_cover` for the calibration and test reads |

## Report levels and cascade

Decided on the calibration split, separately for each mode. The profile settings also decide `catalog_reported`, whether the reads the gene catalog settles for a unit can be reported, on the calibration reads' catalog calls alone.

| Option | Default | What it does |
|---|---|---|
| `profile_min_report_f1` | 0.5 | The `f1_all` a profile unit must reach on the calibration reads |
| `profile_min_eval_positives` | 5 | Distinct positive calibration sequences a profile unit needs before it is judged; below it, `insufficient_evidence` (0: off) |
| `profile_min_cascade_gain` | 0.01 | The smallest calibration gain in `f1_all` that switches the HMM cascade on for a unit |
| `profile_decision_support` | 0.9 | The fraction of bootstrap resamples a level, and the cascade's gain, must hold in |
| `profile_decision_bootstrap` | 1000 | Bootstrap resamples of the calibration source sequences (0: point estimates only; with `profile_min_cascade_gain 0`, the rule up to v10) |
| `annotate_min_report_f1` | 0.5 | The `f1_all` an annotate unit (`--method hmm`) must reach on the calibration proteins |
| `annotate_min_eval_positives` | 5 | Positive calibration proteins an annotate unit needs before it is judged (0: off) |
| `annotate_decision_support` | 0.9 | The fraction of bootstrap resamples an annotate level must hold in |
| `annotate_decision_bootstrap` | 1000 | Bootstrap resamples of the calibration proteins (0: point estimates only) |
| `decision_bootstrap_seed` | 42 | The bootstrap's seed, both modes |

## Gene catalog

| Option | Default | What it does |
|---|---|---|
| `catalog_name` | `gmgc_gut` | The catalog's name in the release (`catalogs/<name>/`, osmotool's `--habitat`) |
| `catalog_source` | GMGC v1.0 human gut, prevalence > 0.5% | Free-text source recorded in `catalog.json` |
| `catalog_source_url` | Zenodo record 10473258 | Source URL recorded in `catalog.json` |
| `catalog_force` | false | Rebuild the catalog annotation even if it exists |

## BacDive and the NaCl model

| Option | Default | What it does |
|---|---|---|
| `bacdive_full` | false | false: the analysis on the 43 pilot strains. true: every BacDive strain with a NaCl maximum and a genome, the salt-in control, and the NaCl model |
| `bacdive_pilot` | the 43 pilot strains | The strain table of the default run |
| `bacdive_snapshot` | none | Harvested BacDive batches (`bacdive_raw/`) to start from, to reproduce a published table |
| `nchunks` | 40 | Tasks per per-genome step (BacDive proteins and annotation) |

## Real-genome validation

| Option | Default | What it does |
|---|---|---|
| `reference_genomes_dir` | none | osmotool's `test_genomes/`, for the system-call check; unset: that check is skipped |
| `test_genomes_dir` | `assets/test_genomes` | Genomes with curated target genes, for the read-recall check |
| `recall_read_pairs` | 20000 | Simulated read pairs per genome for the read-recall check |

## Profiling samples (`-entry profile_samples`)

Only for `-entry profile_samples`, which runs `osmotool profile` on your own samples with a finished release; the build does not use them.

| Option | Default | What it does |
|---|---|---|
| `samples` | none | CSV with columns `id,r1,r2` |
| `db` | none | The release folder to profile with |

## Notes

- **`refs_store`.** Each family's fetch is stored under `<UniProt release>/families/<family>_<key>/`, the key being a hash of the family's entry and `fetch_refs.py`. The fusion-pair candidates and HMM decoys are stored under `<UniProt release>/shared/`. An edited family refetches only itself. A new UniProt release, or any change to `fetch_refs.py`, refetches everything. The cleaning, clustering and split always run.
- **`diamond_min_subject_cover`.** It discards hits that cover only a small part of the reference protein, a failure seen with proX, where reads matched only a fused-in domain. Every short read covers little of a full-length protein, so start low (0.1 to 0.3) and compare `results/metrics/summary.tsv` with a run without it before raising it.
