# Taxonomy split: sensitivity analysis of the v10 benchmark

**Conclusion.** The v10 benchmark (random split) stays the primary result, and
v10's report levels stay as decided. Holding out whole genera costs annotate
about 3 points of precision and profile about 2 points of median f1_all.
Four reported profile units do not hold up on unseen genera (proX, gshF,
betB, proP), and they are named as such rather than re-decided.

Release: `releases/v10-taxsplit` (built and benchmarked 2026-11-16). Log:
`releases/v10-taxsplit_build_benchmark.log`.

## Why

v10 splits each family's positives at random, after CD-HIT clustering at 90%
identity. A close relative of almost every test protein therefore stays in
the training set: the median test protein is 81% identical to its nearest
training reference (DIAMOND blastp, best same-family hit, >= 50% query
cover). Catalog-based tools such as Cayman cluster at 60% before splitting.
The question here is whether v10's numbers depend on that lenient split.

## What changed

Only step 3, `split_train_test.py --split-mode taxonomy`:

- **Positives are split by genus.** The genus is the first word of the
  organism name in the FASTA header (`opuAB|A0A1Y5T4E0|Aquimixticola_soesokkakensis`
  gives *Aquimixticola*). Genera are shuffled and moved whole into the
  held-out set until it reaches the test fraction (20%). At least one genus
  always stays in training. About 20% of genera were held out per family
  (e.g. ectA 128/576, opuBA 12/77).
- **The held-out positives are split into calib and test by genus too**, so
  calibration never sees a genus that the benchmark scores.
- **Hard negatives stay split at random.** They are other genes, so holding
  out genera adds nothing there.

Everything else is the v10 recipe, run from scratch on the new split:

```bash
nextflow run main.nf --release v10-taxsplit --split_mode taxonomy
```

- v10-taxsplit was built on the same fetched, QC'd and clustered sequences as
  v10 (the old pipeline could reuse a release's references), so the split was
  the only difference. The pipeline now fetches its references in every
  build: to compare the two splits again, build both from the same day's
  fetch (`--split_mode random` and `taxonomy`, one after the other).
- New alignments, HMMs and DIAMOND database on the new training set.
- GA cutoffs, DIAMOND cutoffs and the cascade decision recalibrated on the
  new calib split.
- The report levels are v10's, from `families.yaml`. Step 11c warns where
  this split's own calibration would choose differently (below).
- osmotool from branch `report-levels`; the cross-checks agree (3475/3475
  sample x family cells, 3742/3742 sample x unit cells).

## How strict the split is

| | v10 (random) | taxonomy |
|---|---|---|
| Median identity, test protein to nearest training reference | 81.0% | 72.1% |
| Test proteins above 80% | 53% | 24% |
| Test proteins below 60% | 8% | 17% |

Stricter, but not 60%-clustering strict. Orthologs from different genera are
often still quite similar, so a genus holdout does not guarantee low
identity.

## Annotate (full-length proteins, step 11b)

| | v10 HMM | taxonomy HMM | v10 DIAMOND | taxonomy DIAMOND |
|---|---|---|---|---|
| micro precision | 0.921 | 0.890 | 1.000 | 1.000 |
| micro recall | 0.935 | 0.936 | 0.495 | 0.238 |
| micro F1 | 0.928 | 0.912 | 0.662 | 0.385 |
| macro F1 | 0.799 | 0.779 | 0.515 | 0.342 |
| hard negatives called as any family | 7.3% | 10.9% | 2.3% | 2.2% |
| background calls per 10k | 4.0 | 2.0 | 0 | 0 |

HMM annotate recall by identity to the nearest training reference:

| identity | < 40% / no hit | 40-60% | 60-70% | 70-80% | 80-90% |
|---|---|---|---|---|---|
| v10 | 0.83 | 0.76 | 0.91 | 0.95 | 0.97 |
| taxonomy | 0.84 | 0.85 | 0.93 | 0.97 | 0.98 |

- HMM annotate holds up. Recall is unchanged; the cost is about 3 points of
  precision, from more hard negatives being called. Annotate's 0.92/0.94 is
  not an artefact of the random split.
- DIAMOND annotate halves its recall on unseen genera. This supports HMM as
  annotate's default.

## Profile (simulated reads, DIAMOND + cascade, per reporting unit)

f1_all counts every non-background read called as the unit as a false
positive (see `modules/local/benchmark/metrics/templates/compute_metrics.py`).

Over the 27 reported units, median values:

| | v10 | taxonomy |
|---|---|---|
| f1_all, DIAMOND + cascade (what profile reports) | 0.770 | 0.750 |
| precision_all | 0.890 | 0.914 |
| recall | 0.674 | 0.650 |
| f1_all, DIAMOND alone | 0.767 | 0.723 |
| f1_all, HMM alone | 0.673 | 0.649 |

Random-background false positives (DIAMOND + cascade, all families): 2083
per million reads (v10) and 2214 (taxonomy).

Per unit ("test proteins" is the size of the taxonomy test set; "calibration
says" is what this split's own calibration rule would choose, where it
differs from v10's level):

| unit | level | v10 precision_all | v10 recall | v10 f1_all | tax precision_all | tax recall | tax f1_all | Δ f1_all | test proteins (tax) | tax calibration says |
|---|---|---|---|---|---|---|---|---|---|---|
| ectA | gene | 0.998 | 0.674 | 0.805 | 0.986 | 0.612 | 0.755 | -0.050 | 180 | same |
| ectB | gene | 0.821 | 0.616 | 0.704 | 0.705 | 0.65 | 0.677 | -0.027 | 140 | same |
| ectC | gene | 1.0 | 0.671 | 0.803 | 1.0 | 0.634 | 0.776 | -0.027 | 137 | same |
| teaA | gene | 0.879 | 0.669 | 0.76 | 0.898 | 0.644 | 0.75 | -0.010 | 14 | same |
| betL | gene | 0.929 | 0.815 | 0.868 | 0.882 | 0.796 | 0.837 | -0.031 | 79 | same |
| kdpA | gene | 1.0 | 0.955 | 0.977 | 1.0 | 0.94 | 0.969 | -0.008 | 431 | same |
| nhaA | gene | 1.0 | 0.902 | 0.949 | 1.0 | 0.866 | 0.928 | -0.021 | 794 | same |
| proX | gene | 0.605 | 0.504 | 0.55 | 0.382 | 0.378 | 0.38 **< 0.5** | -0.170 | 56 | same |
| proP | gene | 0.374 | 0.727 | 0.494 | 0.309 | 0.638 | 0.417 **< 0.5** | -0.077 | 126 | none |
| otsB | gene | 0.85 | 0.746 | 0.795 | 0.877 | 0.672 | 0.761 | -0.034 | 347 | same |
| mscL | gene | 0.929 | 0.653 | 0.767 | 0.934 | 0.619 | 0.744 | -0.023 | 955 | same |
| galE | gene | 0.82 | 0.899 | 0.858 | 0.819 | 0.881 | 0.849 | -0.009 | 923 | same |
| mazG | gene | 0.982 | 0.81 | 0.888 | 0.979 | 0.803 | 0.883 | -0.005 | 392 | same |
| opu_permease | group | 0.478 | 0.544 | 0.509 | 0.423 | 0.681 | 0.522 | +0.013 | 51 | opuAB gene; opuBB, opuCB none |
| mrpB | gene | 0.979 | 0.564 | 0.716 | 1.0 | 0.506 | 0.672 | -0.044 | 26 | same |
| mrpE | gene | 1.0 | 0.528 | 0.691 | 1.0 | 0.526 | 0.689 | -0.002 | 37 | none |
| mrpF | gene | 1.0 | 0.543 | 0.704 | 1.0 | 0.414 | 0.586 | -0.118 | 52 | same |
| mrpG | gene | 1.0 | 0.515 | 0.68 | 1.0 | 0.472 | 0.641 | -0.039 | 268 | same |
| betA | gene | 0.54 | 0.637 | 0.584 | 0.597 | 0.58 | 0.589 | +0.005 | 168 | same |
| betB | gene | 0.444 | 0.54 | 0.487 | 0.48 | 0.434 | 0.456 **< 0.5** | -0.031 | 153 | none |
| gshA | gene | 0.666 | 0.914 | 0.77 | 0.737 | 0.763 | 0.75 | -0.020 | 171 | same |
| gshF | gene | 0.533 | 0.879 | 0.664 | 0.304 | 0.877 | 0.452 **< 0.5** | -0.212 | 18 | same |
| gshB | gene | 0.89 | 0.906 | 0.898 | 0.914 | 0.892 | 0.903 | +0.005 | 369 | same |
| gor | gene | 0.975 | 0.936 | 0.955 | 0.976 | 0.906 | 0.94 | -0.015 | 116 | same |
| cspA | gene | 0.997 | 0.638 | 0.778 | 0.993 | 0.636 | 0.775 | -0.003 | 194 | same |
| trkA | gene | 0.834 | 0.777 | 0.804 | 0.812 | 0.847 | 0.829 | +0.025 | 255 | same |
| ktrA | gene | 0.812 | 0.745 | 0.778 | 0.98 | 0.672 | 0.797 | +0.019 | 38 | same |

The 15 excluded units stay below 0.5 on the taxonomy split, except otsA
(0.554 -> 0.578) and murB (0.927 -> 0.920); mrpC sits on the bar (0.503 ->
0.499). All three are excluded through `scope: annotate_only`, not through
the report-level rule. Of the units the rule excluded, the closest are mscS
(0.379) and mrpA (0.411).

### Reading the per-unit results

- **Most units barely move:** 22 of 27 change by less than 0.05 in f1_all.
  Precision even rises slightly overall. The loss is mainly recall: reads
  from unseen genera align less often at DIAMOND's thresholds.
- **Four reported units fall below 0.5:**
  - proX (0.55 -> 0.38) and gshF (0.66 -> 0.45) lose precision: more of
    what gets called as them is a paralog or a hard negative. Both drops
    rest on small test sets (56 and 18 proteins), so they are noisy.
  - proP (0.49 -> 0.42) and betB (0.49 -> 0.46) were already just below 0.5
    on the v10 test split.
- **Report-level stability.** This split's own calibration would choose
  differently for proP, betB, mrpE and opu_permease. All were within a few
  hundredths of 0.5 on v10's calibration. The other 36 families keep their
  level. On the taxonomy test reads, opu_permease still reaches 0.522 as a
  group and mrpE 0.689 at gene level.

## Caveats

- Small test sets for several families (teaA 14, gshF 18, mrpB 26, opuBA 12
  test proteins): their per-unit numbers move by chance as well as by split.
- One split, one seed. The size of the drop is an estimate, not a bound.
- Hard negatives are split at random in both releases, so precision against
  hard negatives from unseen genera is not tested here.
- The simulated reads come from back-translated proteins with one fixed
  codon per amino acid. That affects both releases equally, but it means
  nucleotide-level variation between strains is not modelled.

## Decision

- **v10 stays the primary benchmark and release candidate.** Its numbers
  come from a split decided before any test read was scored; this analysis
  is reported next to it.
- **v10's report levels and cascade stay as decided.** Changing them because
  of these test numbers would select on the test split. Deciding levels on
  the taxonomy split's calibration instead would be a legitimate, more
  conservative design, but it is a different design and is not adopted.
- **In the paper:** report annotate under both splits, and name proX, gshF,
  betB and proP as profile units that do not hold up on unseen genera.

## Files

- `releases/v10-taxsplit/results/metrics/summary_by_report_unit.tsv`
  (profile, per reporting unit, all three methods)
- `releases/v10-taxsplit/results/annotate/overall.tsv` and `summary.tsv`
  (annotate)
- `releases/v10-taxsplit/osmo_refdb.profile_report_levels.tsv` (this
  split's own calibration choices, next to v10's levels)
- v10 counterparts under `releases/v10/results/`
