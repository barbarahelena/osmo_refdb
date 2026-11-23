# Decision robustness: the report levels and the cascade under calibration noise

**Conclusion.** The cascade and the report levels are now decided with a
cluster bootstrap over the calibration evidence, not on bare point
estimates. On v10 the cascade's new minimum gain (the gain floor, below)
alone switches three units off: galE, gshA and opu_permease. opu_permease is the one whose calibration
gain reversed on the test reads. In annotate the support test moves one unit,
mrpC, to `none`, and that unit held up on test. The profile bootstrap has
not been run on real reads yet: v10 ships only the calibration tables, not
the reads they were computed from.

Code: `report_rule.py` (the decision helpers), `compute_metrics.py`
(`cluster_bootstrap`, `decide_units`) and `benchmark_annotate.py`
(`protein_bootstrap`, `decide_levels`). Tests: `tests/python/`.

## Why

Both decisions in `decide_units` used to rest on calibration point estimates.
There was no account of how far those estimates move.

- **The cascade** was switched on for any gain in `f1_all`, however small.
  On v10's calibration reads opu_permease gained +0.006 (0.484 to 0.490). On
  the test reads the cascade lowered it, from 0.455 to 0.407.
- **The level threshold** compared `f1_all` with 0.5. A unit whose true value
  sits near 0.5 lands on either side at random. At gene level every member of
  a paralog group must clear it, so one noisy member moves the whole group.

Neither is a peeking problem: the test reads are never used to decide. What
was missing was a minimum effect size and any measure of uncertainty.

## The rule

- **Resampling.** The unit of resampling is the source sequence, not the
  read: reads from one CDS are not independent. Sequences are stratified by
  (family, positive / hard negative), and each stratum is drawn with
  replacement to its own size. Flank reads go with the construct they were
  cut from. The background stays out, as it does from `f1_all`. One draw
  serves every unit and both cascade options, so the on/off comparison is
  paired. Default: 1000 resamples, seed 42.
- **Cascade.** The gain is a unit's calibration `f1_all` with the cascade
  on minus its `f1_all` with the cascade off. The cascade is on at a unit's
  level only if both of these hold:
  - **the gain floor:** the point gain, at 3 decimals, is at least
    `profile_min_cascade_gain` (0.01). v10 had no floor: any gain above 0
    switched the cascade on;
  - **the support test:** the paired gain is positive in at least
    `profile_decision_support` (0.9) of the resamples.
- **Level.** A level needs `f1_all >= profile_min_report_f1` (0.5) both on
  the point estimate and in at least 0.9 of the resamples, at the cascade
  option chosen above. Clearing on the point estimate only gives `none` with
  reason `unstable`. Uncertainty therefore pushes a unit down the order
  (gene, then group, then none). The `profile_min_eval_positives` guard is
  unchanged.
- **Annotate** uses the same level test with settings of its own
  (`annotate_min_report_f1`, `annotate_min_eval_positives`,
  `annotate_decision_support`, `annotate_decision_bootstrap`; the defaults
  are profile's values). Reads and proteins are different benchmarks, so the
  two thresholds are set separately. Each protein is its own source
  sequence, so the resampling is a protein bootstrap. Annotate has no
  cascade.
- **The v10 rule.** `profile_decision_bootstrap = 0` with
  `profile_min_cascade_gain = 0` gives v10's rule. Replayed on v10's calibration tables, it reproduces
  v10's `profile_report_levels.tsv` and `profile_cascade.tsv` cell for cell
  (`tests/python/test_decisions.py`, fixtures in
  `tests/python/data/v10_calibration/`).

## What it changes on v10

### The cascade: the gain floor (measured on v10's calibration tables)

Nine reported units had the cascade on in v10, each because its calibration
gain was above 0. The table applies the gain floor alone to those nine: a
unit keeps the cascade only if its calibration gain is at least 0.01. The
support test is not in the table, because it needs the calibration reads
(see below). It can only switch more cascades off.

- **calibration gain:** `f1_all` cascade on minus off on the calibration
  reads, the number the decision uses.
- **test gain:** the same difference on v10's test reads
  (`summary_by_report_unit.tsv`, `diamond_cascade` against `diamond`).
  It is shown for reading only and was never used to decide.
- **v10:** the cascade under v10's rule (on for any gain above 0).
- **floor 0.01:** the cascade under the gain floor alone.

| unit | calibration gain | test gain | v10 | floor 0.01 |
|---|---|---|---|---|
| ectB | +0.039 | +0.004 | on | on |
| teaA | +0.280 | +0.232 | on | on |
| proX | +0.336 | +0.317 | on | on |
| galE | +0.001 | +0.002 | on | **off** |
| opu_permease | +0.006 | −0.048 | on | **off** |
| betA | +0.155 | +0.140 | on | on |
| betB | +0.027 | +0.004 | on | on |
| gshA | +0.002 | +0.002 | on | **off** |
| trkA | +0.109 | +0.104 | on | on |

Of the nine cascaded units, only opu_permease lost `f1_all` on test. The
other two that the floor switches off gained +0.002 on test, so switching
them off costs that little and saves their hmmsearch time. Two gains just
above the floor, ectB +0.039 and betB +0.027, shrank to +0.004 on test.
Whether the paired support keeps them on is one of the profile numbers still
to measure (below).

### Annotate: the level test (measured on v10's calibration proteins)

The run used osmotool 0.10.0 annotate (`--method hmm`) on positive.calib
against negative.train: 38,558 proteins, 1000 resamples.

- One unit changes. mrpC has `f1_all` 0.548 on 48 calibration proteins, with
  p10 0.493 and p90 0.602, and 0.886 of resamples clear 0.5. It becomes
  `none` for `unstable`.
- Every other unit keeps its v10-recipe level: 29 gene, opu_permease as a
  group, and 12 none, mrpC included.
- On the test split mrpC scores `f1_all` 0.629 (39/48 positives found). This
  is the cost the rule accepts by design: when it is unsure, it reports
  nothing.

### Profile: the level test and the cascade support (not yet measured)

Both need the calibration reads, which are not in the release. Two facts
bound what they can change:

- `families.yaml` sets `report_level` for every family, so in profile the
  level test cannot change what is reported. It changes `calibration_level`,
  and a unit families.yaml keeps where the rule now says `none` becomes an
  `override`. The units nearest 0.5 on v10's calibration reads are:
  - proP 0.519 and proX 0.527, both gene;
  - mrpC 0.476 (`annotate_only`);
  - two that are already overrides: betB 0.490 and opu_permease 0.490.

  proX and proP are also two of the four units that did not hold up on unseen
  genera (`taxonomy_split.md`).
- The support test can only switch cascades off: of the six left on after
  the floor, those whose paired gain is positive in fewer than 0.9 of the
  resamples.

`results/calib_metrics/decision_sensitivity.tsv` lists, per unit, the
decision on point estimates alone next to the rule's. The next build writes
it.

## Cost

These timings are from this machine. The cluster run will differ.

- **Annotate:** 1000 resamples of 38,558 proteins take 7.2 s in pure Python
  (Python 3.9). The scorer's container keeps no numpy.
- **Profile:** a synthetic calibration table at v10's size takes 5.1 s for
  1000 resamples, with 0.85 GiB peak RSS including the table. The table had
  0.8M reads over 38,558 source sequences, 51 units, and both cascade options.
  The time grows with the number of sequences each unit's false positives
  come from. It is bounded by the number of sequences, but the real
  distribution has not been measured.
