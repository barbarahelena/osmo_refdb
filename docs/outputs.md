# Outputs

Everything is published into the release folder, `--outdir`. The table lists each process that publishes something, where it goes, and what it is. Processes not listed pass their results on without publishing them. Where they publish is set in `conf/modules.config`.

## The database (what osmotool reads)

| Process | Published to | What it is |
|---|---|---|
| `DIAMOND_MAKEDB` | `osmo_refdb.dmnd` | The DIAMOND database |
| `DIAMOND_TRAINREFS` | `osmo_refdb.train_refs.faa` | The sequences the DIAMOND database is built from |
| `DIAMOND_CALIBRATE` | `osmo_refdb.diamond_cutoffs.tsv` | Per-family minimum DIAMOND bitscores |
| `DIAMOND_SCOPE` | `osmo_refdb.profile_excluded_families.txt`, `osmo_refdb.annotate_excluded_families.txt` | Families and decoys kept out of reported output |
| `HMM_CALIBRATE` | `hmms/` | One HMM per family and decoy, the pressed `osmo_refdb.hmm`, the GA cutoffs (`cutoff_manifest.tsv`) and the calibration scores |
| `HMM_CALIBRATE` | `qc_scorecard.tsv` | Per-family QC flags: negative purity, length outliers, HMM cutoff overlap |
| `BENCHMARK_METRICS_CALIB` | `osmo_refdb.profile_cascade.tsv` | Which units get the HMM cascade in `osmotool profile`, and at which threshold |
| `BENCHMARK_METRICS_CALIB` | `osmo_refdb.profile_report_levels.tsv` | The level each unit is reported at in `osmotool profile`, with the calibration numbers behind it, and whether the reads the gene catalog settles for it can be reported (`catalog_reported`; osmotool does not read it yet) |
| `BENCHMARK_ANNOTATESCORE_CALIB` | `osmo_refdb.annotate_report_levels.tsv` | The same for `osmotool annotate --method hmm` (osmotool does not read it yet) |
| `MARKERS_PACKAGE` | `markers/` | The marker genes' DIAMOND database, its source FASTA, excluded labels and manifest |
| `CATALOG_PACKAGE` | `catalogs/<name>/` | The gene catalog: annotation, bwa-mem3 index, hmmsearch domain table, `catalog.json`, `SHA256SUMS` |
| `CATALOG_SCREEN` | `catalogs/<name>/screen.tsv` | Families that claim catalog genes that look like something else |
| `NACL_MODEL` | `nacl_model/` | The NaCl-max model for osmotool (`--bacdive_full` only) |

## References and alignments

| Process | Published to | What it is |
|---|---|---|
| `SPLIT_REFS` | `refs/` | Per family: positives and negatives, clustered and split (`positive.train/calib/test.faa`, `negative.train/test.faa`), and the reference manifests (purity, length outliers, RefSeq, random background) |
| `DIAMOND_EXTRAREFS` | `refs/` | Decoy (`<family>.decoy.faa`) and fused-ORF (`<a>_<b>_fused.faa`) references |
| `HMM_DECOYSELECT`, `CDHIT_DECOY` | `refs/` | The HMM decoys' input and clustered sets |
| `MAFFT_ALIGN`, `TRIMAL` (and their decoy aliases) | `alignments/` | Per family and decoy: the alignment (`.aln.fasta`) and the trimmed alignment (`.trimmed.fasta`) |

## Benchmark results

| Process | Published to | What it is |
|---|---|---|
| `BENCHMARK_METRICS_CALIB` | `results/calib_metrics/` | Calibration reads: `summary.tsv`, `summary_by_group.tsv` (DIAMOND, HMM, the cascade, and the catalog calls alone as `catalog_only`), `best_threshold_summary.tsv`, `decision_sensitivity.tsv`, ROC/PR plots |
| `BENCHMARK_METRICS_TEST` | `results/metrics/` | Test reads: `summary.tsv`, `summary_by_group.tsv`, `summary_by_report_unit.tsv` (the numbers for the paper), `average_precision.tsv`, `verification.tsv`, `best_threshold_summary.tsv`, ROC/PR plots; with the options, `summary_by_read_length.tsv` and `background_fpr*.tsv` |
| `BENCHMARK_ANNOTATESCORE_CALIB` | `results/annotate_calib/` | Calibration proteins: osmotool annotate's output and calls per method, `summary.tsv`, `summary_by_group.tsv` |
| `BENCHMARK_ANNOTATESCORE_TEST` | `results/annotate/` | Test proteins: the same, plus `average_precision.tsv` and `verification.tsv` |
| `BENCHMARK_CATALOGSCORE` | `results/catalog_benchmark/` | The catalog benchmark: catalogs A and B scored per reporting unit |

## Validation

| Process | Published to | What it is |
|---|---|---|
| `BACDIVE_ANNOTATE` | `phenotype_validation/annotate/`, `phenotype_validation/strategy_control/annotate/` | osmotool annotate per BacDive genome |
| `BACDIVE_MATRIX` | `phenotype_validation/feature_matrix.tsv` | Strain x family gene counts with each strain's phenotype |
| `BACDIVE_ANALYZE` | `phenotype_validation/` | The phenotype analysis: `family_ranking.tsv`, `null_panels.tsv`, `cv_results.tsv`, the salt-in control, and each test's report (`*.txt`) |
| `NACL_MODEL` | `nacl_model_build/` | The model's training table, cross-validation and SHAP explanation (`--bacdive_full` only) |
| `VALIDATION_CHECKSYSTEMS` | `real_genome_validation/reference_systems_check.*` | System calls on the reference genomes against the expected systems |
| `OSMOTOOL_PROFILE_RECALL`, `OSMOTOOL_PROFILE_RECALL_CATALOG`, `HMMER_HMMSCAN` | `real_genome_validation/` | Per genome: osmotool profile's counts and alignments, without and with the catalog, and hmmscan on the reads |
| `VALIDATION_READRECALL` | `real_genome_validation/` | Per genome: the read truth and the recall/precision comparison |

## BacDive data (release-independent)

| Process | Published to | What it is |
|---|---|---|
| `BACDIVE_HARVEST` | `<data_dir>/` | The BacDive strain table, classes, the salt-in control's strains, and the harvested batches (`bacdive_raw/`); `--bacdive_full` only |
| `POOL_PFAM`, `POOL_EXTRACT` | `<data_dir>/random_panel/` | The random-panel pool's models (`pool.hmm`) and the Pfam release they come from |

## Provenance

| Process | Published to | What it is |
|---|---|---|
| `PROVENANCE` | `PROVENANCE.json` | Code commit, containers and tool versions, inputs, parameters, checksums of every output, licence |
| `PROVENANCE` | `software_versions.yml` | Tool versions per process |

Nextflow's own reports (trace, timeline, run report) go to `pipeline_info/` in the launch directory. `-entry profile_samples` publishes each sample's osmotool profile output into `--outdir`.
