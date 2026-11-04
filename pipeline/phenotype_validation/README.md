# Phenotype validation

The `benchmark` step and `real_genome_validation/` both ask the same kind of
question: can DIAMOND/HMM **find** family X. Neither asks whether the gene set
the panel recovers means anything physiologically. This directory is the third
tier, and it tests a different claim:

> does a genome's osmo_refdb profile carry information about that organism's
> measured salt tolerance, beyond what genome size and phylogeny already tell
> you?

The endpoint is an externally measured phenotype — growth at a stated NaCl
concentration, curated by DSMZ — not a proxy derived from the same annotation
layer the database was built from. That independence is the whole point. A
breadth score computed from our own annotations and scored against a label
derived from those same annotations measures internal consistency and cannot
fail.

This workflow does **not** replace the read-level benchmark. That one validates
detection and calibrates cutoffs. This one validates the panel's biological
content. Different jobs.

## Status

Pilot complete on a 43-strain random sample. Full BacDive ID sweep harvested
and the phenotype table built (2026-11-04); genome download, annotation and
analysis not yet run. The
analysis plan below is **pre-registered** — written before the full harvest, so
that the salt-in prediction and the blocking strategy are commitments rather
than post-hoc descriptions. Amendments belong in separate commits with a stated
reason.

## Prerequisite: osmotool pin

**This workflow needs osmotool at `2826ebc` or later.** The pilot was run
against that commit (verified from the image's `direct_url.json`, not assumed).
Earlier pins, including the `4cbc3ba` on `main` at the time this branch was cut,
predate `quantifier.KNOWN_FAMILIES` being corrected from an original ~15-family
list to the full 43-family panel — so families like `mrpG`, `gshB` and `trkA`
are not in the old pin's known-family list at all, and `annotate` would silently
under-report.

`main` now pins osmotool `v0.7.0` (PR #15), which satisfies this. v0.7.0 also
switches `annotate`'s HMM step from hmmscan to hmmsearch (osmotool 2358bc0),
so the pilot's throughput figures, measured on `2826ebc` with hmmscan, are
conservative for the full run. Rebuild the image before running, and check:

```bash
grep 'osmotool.git@' Dockerfile    # expect v0.7.0 or later
```

## Data source: BacDive

`Physiology and metabolism > halophily` in the BacDive v2 API:

```json
{"salt":"NaCl","growth":"positive","tested relation":"growth",  "concentration":"3-25 %"}
{"salt":"NaCl","growth":"positive","tested relation":"optimum", "concentration":"10 %"}
```

`tested relation` takes growth / optimum / maximum / minimum, so a tolerance
interval and an optimum are both available. Negatives are encoded
`growth: "no"` — note the token, it is not `"negative"` — and give an upper
bound on true tolerance, which makes the endpoint interval-censored rather
than a point value.

Registration was dropped in February 2026; `/v2/fetch/{id}` is open and accepts
100 semicolon-separated IDs per call. Use **v2**: v1 still responds but is
frozen at April 2025 data.

Coverage, measured on a seeded random sample of 1000 IDs (558 resolved):

| | share of resolved |
|---|---|
| total strains in BacDive | ~97,650 projected (matches DSMZ's published count) |
| halophily record | 10.9% |
| genome accession | 20.8% |
| both | 8.2% |
| genome + parseable quantitative NaCl max | 7.0% → **~6,800 strains** |

Full ID sweep (2026-11-04, IDs 1–182,000 in 1,820 batches), which replaces the
sample's projections with counts:

| | count |
|---|---|
| strains in BacDive | **102,187** (IDs end at ~179,200) |
| halophily record + genome accession | 6,643 |
| … with a quantitative NaCl max | **6,344** |
| … with a no-growth upper bound | 1,794 |
| distinct genera (quantitative max) | 1,608 (3.95 strains/genus) |
| archaea | 110 |
| pre-registered salt-in strains | 59, in 23 genera |
| assembly level (quantitative max) | 1,024 complete, 115 chromosome, 2,718 scaffold, 2,487 contig |

The sample's ~97,650 total undercounted because it drew from IDs 1–175,000,
the same cap the sweep originally had. 163 NaCl records are skipped by
`build_table.py` rather than parsed: 99 in molar units and 64 in g/L.

Strain-level matters: tolerance varies within a species, and records carry
`INSDC accession` plus `assembly level` and a strain↔genome match `score`, so
the genome is matched to the strain that was actually assayed. A species-level
join cannot do this.

## Why not the Madin trait compilation

`condensed_traits_NCBI.csv` was evaluated and rejected. Its only salinity field,
`range_salinity`, is categorical rather than a concentration (the numeric
`d1_*`/`d2_*` columns are cell dimensions), covers 684 species, and derives
entirely from PATRIC's curator metadata field — one raw entry is literally
`4%`. Decisively, their `data/conversion_tables/renaming_range_salinity.csv`
maps `Mesophilic` → `moderate-halophilic`: a temperature trait folded into the
dominant salinity class, and at 268 raw entries against 167 genuinely moderate
ones it accounts for most of that class. Madin's BacDive input carries only
`metabolism`, so routing through it yields none of BacDive's NaCl data.

Re-deriving from `output/prepared_data/patric.csv` and dropping `Mesophilic`
salvages 494 species on 3 effective classes. BacDive's ~6,800 strains on a
continuous endpoint is the better instrument by a wide margin.

## Pre-registered analysis plan

### Baselines

The panel is never judged on its own correlation. It is always compared against:

- `total_proteins` (proteome size, from the `gene_counts.tsv` header)
- genome size and GC
- a size-matched random panel of 43 non-osmotic families

The defensible claim is "the panel beats these", never a bare correlation or AUC.

### Two analysis sets

**Set A — unbiased.** Every strain with a genome and a quantitative NaCl max,
from the full ID sweep. Headline estimate, evaluated under genus-grouped
cross-validation (all strains of a genus in one fold) and again blocked at
family level. Report both; the gap between them is informative.

**Set B — ancestry-controlled.** Genera with ≥3 strains spanning ≥4 percentage
points of NaCl max. Grouped CV alone does not fix confounding — a clean
across-clade correlation can still be entirely ancestry. Only within-genus
contrast addresses it. Model: `nacl_max ~ panel + (1|genus)`, so the fixed
effect is estimated from within-genus variation.

How many genera qualify is unknown until the sweep, which is why the sweep runs
*before* the stratification is chosen. The pilot's random sample gave ~1 strain
per genus (37 genera / 41 strains), so Set B requires deliberate targeting and
cannot be assembled from a random draw.

**Amendment (2026-11-04).** Written after the sweep and the phenotype table,
before any panel output on the sweep's strains existed. The full sweep is not a
random draw, and it already holds 324 genera with ≥3 strains spanning ≥4
percentage points of NaCl max (4,291 strains; largest *Streptomyces* 337,
*Pseudomonas* 148, *Flavobacterium* 119, *Paenibacillus* 118, *Vibrio* 86). Set B
is therefore built from the sweep as-is, with no separate targeting step. The
Set B criterion itself is unchanged.

Quantify the confound directly as well: ANOVA of `nacl_max` by
genus/family/order gives the fraction of phenotype variance explained by
taxonomy alone. If that is most of it, report it.

**All three ancestry controls need multi-strain genera to mean anything**, which
makes this a harvest requirement and not just an analysis choice. At ~1 strain
per genus, eta² is forced toward 1.0 by construction (between-group variance
*is* total variance), and the genus-cluster bootstrap degenerates into the
strain bootstrap because resampling clusters of size 1 is resampling strains.
The pilot hits exactly this: eta² by genus = 0.96 at 1.1 strains/genus, which is
arithmetic, not phylogenetic signal. `analyze.py` detects the condition and
labels those outputs DEGENERATE rather than printing them straight.

**Diagnostic, committed in advance:** if within-genus rho stays positive and the
top-ranked families remain K⁺/Na⁺/compatible-solute transporters, the signal is
real. If within-genus rho collapses toward zero while across-genus rho stays
high, it was ancestry.

### Salt-in stratum

Membership is fixed a priori in `salt_in_taxa.tsv`, assigned from published
physiology before any panel output was inspected. Salt-in organisms accumulate
molar KCl rather than organic osmolytes, so they reach maximum NaCl tolerance
with the panel's primary mechanism absent. That makes them a **negative
control**, not merely a confounder: if the panel fires on them, it is a false
positive that no other stratum would reveal.

The prediction is family-specific, and deliberately not "low breadth" — bulk K⁺
accumulation *is* the salt-in mechanism, so the K⁺ machinery should be present:

| family group | prediction in salt-in strains |
|---|---|
| K⁺ uptake (`ktrA`, `trkA`) | present / high |
| compatible solute synthesis + transport (`ectABC`, `betL`, `opuA`, `otsAB`, `proP`) | low / absent |

Treatment, decided up front: primary analysis on salt-out strains only, with
salt-in reported separately as a predicted-divergence cohort; secondary
analysis over everything with a strategy indicator and interaction term.
Primary-on-salt-out keeps the monotone gene-dose assumption honest; dropping
them without declaring it would be cherry-picking.

Haloarchaea are nearly absent from a random draw (1 of 46 archaeal strains in
the pilot), so this stratum needs a targeted `/v2/taxon/{genus}` pass over the
genera in `salt_in_taxa.tsv`.

**Amendment (2026-11-04).** Same timing as the Set B amendment above. The
targeted pass is dropped as redundant: the full sweep fetched every BacDive
ID, so a taxon query can only return strains already harvested. The sweep
gives 59 salt-in strains across 23 of the genera in `salt_in_taxa.tsv`.
Membership and the family-specific prediction are unchanged.

### Power

Sample size is set from a **pre-specified** effect size — the smallest rho worth
acting on — never from an observed pilot correlation. Feeding a pilot's own
estimate into a power calculation is circular and biased: a small sample that
looks impressive is usually an overestimate, so it flatters itself and
recommends stopping early. For reference, n ≈ 85 for rho = 0.3, n ≈ 196 for
rho = 0.2, n ≈ 785 for rho = 0.1.

## Pilot result (n = 43 strains, 41 with a quantitative max)

Mechanics, all clean: 43/43 accessions resolved to a downloadable assembly with
no suppressed records, 43/43 annotated with `rc=0`, zero unparsed concentration
strings. Throughput 7.28 s/genome at 10 concurrent single-threaded jobs,
projecting ~13.7 h and ~35 GB for 6,800 strains.

Smoke test only — at n=41 with 43 features these numbers carry no weight, and
exist to confirm the chain works and the signs are not absurd:

```
breadth (families present)      rho=+0.375  [+0.065, +0.637]
total copies                    rho=+0.312  [+0.006, +0.576]
BASELINE total_proteins         rho=-0.242  [-0.507, +0.070]

breadth | total_proteins        rho=+0.508  [+0.222, +0.713]
copies  | total_proteins        rho=+0.553  [+0.293, +0.732]
```

Genome size does not predict the phenotype here and is mildly negative, so it
acts as a suppressor — which is why partialling it out raises the panel's
correlation rather than lowering it.

More interesting than the correlations is the family ranking: `ktrA`, `mrpD`,
`betL`, `mrpG`, `mrpF` at the top, `kdpA` most negative. `kdpA` is the
high-affinity K⁺ scavenger for K⁺-poor environments, whereas bulk K⁺ in
halotolerant organisms moves through lower-affinity, higher-flux Trk/Ktr
systems — so `ktrA` up and `kdpA` down is two K⁺ strategies separating in the
mechanistically correct direction.

Three deflations, weighted heavily:

1. `mrpD`/`mrpF`/`mrpG` are subunits of one antiporter complex and co-occur by
   construction. Three of the top five is **one** observation, not three.
2. 37 genera across 41 strains: the phylogenetic confound is entirely
   unaddressed and the whole pattern could be clade structure.
3. Low-breadth strains track proteome size (*Campylobacter*, 1623 proteins → 3
   families). Since genome size runs against the phenotype here this works
   against the result, not for it.

## Workflow

Scripts live here; output goes to `releases/<release>/phenotype_validation/`
(git-ignored), following the `real_genome_validation/` convention.

```bash
PV=releases/v8/phenotype_validation_full   # pilot: releases/v8/phenotype_validation

# 1. harvest BacDive: `sweep` (every ID, ~1,820 calls) or `sample` (the
#    pilot's seeded 1000 IDs). v2 API, 100 IDs per call, no credentials.
#    Keep sample and sweep in separate dirs; the script refuses to mix them.
bash  pipeline/phenotype_validation/harvest_bacdive.sh "$PV" sweep

# 2. phenotype table: one row per strain with NaCl + a genome accession
python pipeline/phenotype_validation/build_table.py --dir "$PV"

# 3. genomes, resolving unversioned GCA accessions via NCBI's FTP listing
bash  pipeline/phenotype_validation/fetch_genomes.sh "$PV"

# 4. real osmotool annotate in the pinned image (10x1 is the fast config)
docker run --rm -v "$PWD":/work -w /work --entrypoint bash osmo_refdb:latest \
    pipeline/phenotype_validation/run_annotate.sh "$PV" releases/v8 10 1

# 5. join features to phenotype, then read out
python pipeline/phenotype_validation/build_matrix.py --dir "$PV"
python pipeline/phenotype_validation/analyze.py      --dir "$PV"
```

## Notes and gotchas

- **`--threads` on `annotate` is counterproductive.** Measured on a
  5330-protein genome: hmmscan takes 13 s at `--cpu 1` and 16 s at `--cpu 8`.
  Threading overhead exceeds the gain. Use 1 thread per job and parallelise
  across genomes instead. (Measured with hmmscan on `2826ebc`; not re-measured
  for hmmsearch under v0.7.0.)
- **hmmsearch is 4.4× faster than hmmscan with an identical hit set**, verified
  across all 43 pilot genomes (526 s → 119 s, 43/43 agree, 0 differ). This is
  safe *specifically because* osmotool uses `--cut_ga`, a bit-score gathering
  threshold; E-values differ between the tools, since hmmscan's search space is
  43 models while hmmsearch's is thousands of sequences. If the cutoff were ever
  changed to an E-value, this equivalence breaks. osmotool v0.7.0 adopted
  hmmsearch (2358bc0); this comparison was a manual benchmark run before that.
- With hmmsearch in osmotool v0.7.0, Prodigal (~11 s, single-threaded) should
  dominate `annotate` and set the per-genome floor.
- Effective parallel speedup is only ~3.2× across 10 cores (~32% efficiency) on
  an M4, whose cores are 4 performance + 6 efficiency. More concurrency will not
  help much; the remaining win is algorithmic.
- Memory is irrelevant: a running `annotate` job uses ~16 MiB.
- The feature matrix has 44 columns, one more than `families.yaml` has entries.
  The extra is `mrpA_mrpB_fused`, from `08d_build_fusion_refs.py`. Not a bug.
- Genus *names* are a coarse proxy for phylogenetic distance and NCBI/LPSN
  genera are unevenly split. A stricter version of the blocking in Set A/B uses
  16S identity or tree distance; BacDive records carry 16S accessions, so the
  data for that is already available.
- Do not re-run the full ID sweep casually. It is ~1,820 calls against a public
  academic service; cache `bacdive_raw/` and reuse it. The 2026-11-04 sweep
  took ~1.5 h at ~3.4 s/batch; 7 batches failed (proxy 502s, 180 s timeouts,
  HTTP/2 framing errors) and a second run fetched them, since cached batches
  are skipped.
