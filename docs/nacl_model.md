# NaCl-max model

This model predicts the highest NaCl concentration (% w/v) at which a salt-out prokaryote grows,
from `osmotool annotate` output alone. It is trained here on the BacDive phenotype-validation set
(`modules/local/bacdive/analyze/templates/`). It ships in the release as `nacl_model/` so that osmotool can
apply it to new genomes.

The model is tied to the release it was trained on: family counts depend on that release's HMM
cutoffs. The current model was trained on **v10** and must be used with
`osmotool annotate --method hmm` against v10.

## Files

| file                         | role                                                                                                                                                                      |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `features.py`              | Builds the features from `<prefix>.gene_counts.tsv` / `.systems.tsv`. Stdlib only; osmotool's reference.                                                              |
| `train.py`                 | Training table, hyperparameter search, nested CV, ablations, final model, intervals, salt-in detector, `feature_schema.json`.                                           |
| `explain.py`               | SHAP importance, stability, interactions, ridge reference, plots, and the README shipped with the model.                                                                  |
| `predict.py`               | Reference inference on annotate output (osmotool's test oracle).                                                                                                          |
| `../package_nacl_model.py` | Copies the build into `<release>/nacl_model/` with checksums and `self_check.json`: the model's raw output for 5 training genomes, which osmotool recomputes on load. |

## Running

```bash
# environment: xgboost >= 3.1, scikit-learn, pandas, scipy, matplotlib (all in environment.yml).
# xgboost < 3.1 loads the saved model without error but misreads its base score
# (0.48 instead of 5.92), so every prediction is off; osmotool and predict.py refuse it.
python3.11 -m venv .venv-nacl && .venv-nacl/bin/pip install xgboost scikit-learn pandas scipy matplotlib

# in the pipeline: NACL_MODEL, with --bacdive_full
nextflow run main.nf --outdir <v10 folder> --release v10 -resume --bacdive_full
# = train.py (~20 min) + explain.py (~10 min) + package_nacl_model.py

# predict, on osmotool >= 0.8.0 annotate output (header has gc and n_contigs)
.venv-nacl/bin/python modules/local/nacl/model/templates/nacl_model/predict.py <v10 folder>/nacl_model out/genomeA out/genomeB
```

On macOS, set `MPLCONFIGDIR` to a writable directory if `~/.matplotlib` is not writable.

## Design

**Features** (46). Everything comes from the annotate output, with no taxonomy:

- `total_proteins`, `gc`, `n_contigs`;
- the `raw_count` of the 43 osmo families that occur in the training set (`mrpA_mrpB_fused` and
  `opuBC` are zero in every Set A strain).

Tested and left out, because they add nothing in genus- or family-blocked CV (`cv_results.tsv`):

- `genome_mb`: also, osmotool counts N's in `genome_bp` and the training data does not;
- the breadth summaries `n_families_present` and `total_osmo_copies`;
- the `systems.tsv` calls.

Leaving them out also stops correlated features from splitting the SHAP credit.

**Target.** `nacl_max_pct` of the 6,264 BacDive salt-out strains with a quantitative value. Salt-in
strains are not in the training set; they are covered by a separate warning (see below).

Tested and not adopted:

- a sqrt target: higher Spearman, higher RMSE;
- 1/sqrt(n_genus) weights;
- dropping the 326 censored labels.

None of these changed the results meaningfully.

**Algorithm.** XGBoost regression trees; the hyperparameters come from a random search scored by
genus-blocked CV. XGBoost fits about as well as sklearn's HistGradientBoosting (the existing CV in
`analyze_extended.py`), and clearly better than ridge. It was chosen for two practical reasons:

- it computes exact TreeSHAP itself (`pred_contribs`), so osmotool needs no `shap` package;
- it stores the model as portable JSON instead of a pickle.

**Validation.** Nested CV:

- outer loop: 5 folds blocked by genus, and separately by taxonomic family;
- inner loop: the hyperparameter search is rerun within each training fold.

The numbers therefore apply to genomes from genera (or families) the model has never seen.

**Uncertainty.** The intervals are split-conformal: the 10th/90th percentiles of the
genus-blocked OOF residuals, in 4 bins of the predicted value.

**Salt-in warning.** A logistic regression (salt-in vs salt-out, 84 vs 6,550 strains, same
features) sets `flags=salt_in_like` when p ≥ 0.5.

## Results (v10, 2026-11-15)

|                                                                              | genus-blocked            | family-blocked |
| ---------------------------------------------------------------------------- | ------------------------ | -------------- |
| Spearman                                                                     | 0.686                    | 0.641          |
| MAE (% NaCl)                                                                 | 2.52                     | 2.73           |
| R²                                                                          | 0.49                     | 0.39           |
| HGB baseline (Spearman)                                                      | 0.681                    | 0.635          |
| ridge baseline (Spearman)                                                    | 0.639                    | 0.620          |
| baseline model, total_proteins + gc only, own nested search (Spearman / R²) | 0.228 / 0.06             | 0.111 / −0.03 |
| 80% interval coverage                                                        | 80.0% (calibrated on it) | 76.3%          |

- **Range compression.** Per bin of the true value (genus-blocked):

  - strains with ≤1% average a prediction of 3.1%;
  - strains with >15% (mean 21.9%) average 13.2%.

  So the model works as a ranking plus a rough level. It does not report a tolerance class,
  because one derived from the prediction would be wrong for most low and extreme strains.
- **What it weighs** (mean |SHAP|, in % NaCl):

  - `betL` (0.85; more copies → higher);
  - `total_proteins` (0.58; larger genome → lower);
  - `kdpA` (0.48; present → lower);
  - `ectA` / `ectC` (0.40 / 0.26; higher);
  - the K+ uptake systems trkH, ktrA and trkA (higher);
  - the mrp subunits (higher).

  Per group, the strongest is compatible-solute uptake, then synthesis, K+ uptake, genome
  size/GC and Na+/H+ antiport. The ranking is stable over the genus folds (ρ 0.92–0.97). It
  matches the within-genus correlations of the BacDive validation (betL, ectABC, mrp, kdpA
  negative). See `<release>/nacl_model/README.md` and `explain/`.
- **Salt-in strains** (72 with a label): median true 29.6%, median prediction 4.9%. This is out
  of domain, and the detector reaches precision 0.87 and recall 0.96 in genus-blocked CV. Its 12
  false alarms are:

  - 4 non-halophilic archaea (Pyrococcus, Palaeococcus, Staphylothermus, Archaeoglobus);
  - anaerobic Firmicutes (Sporosalibacterium, Fusibacter, Alkaliphilus, Natranaerobius);
  - Rubrobacter (2);
  - Gramella;
  - Longibacter, a Salinibacter relative.
- **Parity.** osmotool 0.8.0 on raw genomes (20 BacDive strains, `releases/v10/nacl_model_build/parity/`)
  gives features identical to the training table: protein count, GC, contig count and all family
  and system counts.
- **Reproducibility.** Refitting with the same seed gives a byte-identical model.
- **Test genomes** (osmotool `test_genomes`; their genera are in the training data):

  | genome           | prediction               |
  | ---------------- | ------------------------ |
  | Collinsella      | 1.8%                     |
  | E. faecalis      | 5.6%                     |
  | E. coli          | 7.2%                     |
  | B. subtilis      | 8.2%                     |
  | L. monocytogenes | 10.7%                    |
  | S. aureus        | 11.1%                    |
  | H. elongata      | 13.3%                    |
  | H. volcanii      | flagged `salt_in_like` |
