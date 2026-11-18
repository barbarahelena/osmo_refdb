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
MAFFT, cd-hit, InSilicoSeq, etc. installed just to run `osmotool`.

`osmotool profile` (short reads) uses DIAMOND as the primary caller.
Since osmotool v0.7.0 it also runs a DIAMOND+HMM cascade by default:
reads DIAMOND assigns to a family flagged in
`osmo_refdb.profile_cascade.tsv` (families with low DIAMOND benchmark
precision) get a second opinion from HMMER, using a short-read bitscore
threshold rather than the full-length GA cutoff. So `profile` needs
HMMER too unless you pass `--no_cascade`. Since v10 a release also ships
`osmo_refdb.profile_report_levels.tsv`, the level at which `profile`
reports each family: on its own (`gene`), summed with the close paralogs
short reads can't tell it from (`group`), or not at all (`none`, listed as
excluded rather than left out). Both files are decided together on
calibration reads, never on the benchmark reads (see steps 9c-11c). `osmotool annotate` (full-length ORFs from
an assembly) defaults to HMM (`--method hmm`) since real-genome testing
showed it meaningfully more specific than DIAMOND for this reference
database (see `hmms/osmo_refdb.hmm` + `--hmm_db`); pass `--method diamond`
if you want to avoid the HMMER dependency there too.

## Quick start

The build is a Nextflow workflow (`main.nf`; see `docs/hpc.md` for the cluster setup and
each module's `container` line and `environment.yml` for its image):

```bash
nextflow run main.nf --release v1                           # build, benchmark, catalog, screen; stops at the gate
nextflow run main.nf --release v1 -resume --gate_passed true   # BacDive, NaCl model, validation
```

## What the workflow does

`RELEASE_BUILD` (`subworkflows/local/release_build`), one process per step, each with its
scripts in `modules/local/<step>/templates/`: fetch refs, split, align, hmmbuild + calibrate,
DIAMOND database, then the calibration and test benchmarks. `RELEASE_CATALOG` adds the gene
catalog and the screen; after the gate come the BacDive analysis, the NaCl model and the
real-genome validation. `families.yaml` is the input (`--families_file` to point at a subset).

Everything for a given release name is written to `releases/<name>/`
(git-ignored), so you can keep multiple versions around and reproduce any
of them independently.

**Not yet converted from the old `run_pipeline.sh`:** the taxonomy split and read-length /
background options of the benchmark below. The sections that follow still describe them with the old
`run_pipeline.sh` commands.

By default, step 3's train/test split is a random shuffle of sequences —
a held-out test sequence is still likely to have a close relative in
train. To instead hold out whole genera (testing detection of genuinely
divergent homologs, the scenario where HMMs classically have an edge over
pairwise search), run with `SPLIT_MODE=taxonomy`:

```bash
SPLIT_MODE=taxonomy bash run_pipeline.sh build v2-taxonomy-split
```

The benchmark step also accepts:

```bash
READ_LENGTHS="100,150,250,300" \
BACKGROUND_R1=/path/to/real_metagenome_R1.fastq.gz \
BACKGROUND_R2=/path/to/real_metagenome_R2.fastq.gz \
bash run_pipeline.sh benchmark v1
```

`READ_LENGTHS` simulates at each exact length (forces wgsim for exact
control) and reports precision/recall/F1 broken out per length bucket in
`results/metrics/summary_by_read_length.tsv`. `BACKGROUND_R1`/`BACKGROUND_R2`
run a real paired-end metagenome sample with no known osmoadaptation
content through the same DIAMOND/HMM benchmark and report its
false-positive rate in `results/metrics/background_fpr.tsv` — a more
realistic FPR estimate than UniProt-derived synthetic hard negatives alone.

`DIAMOND_MIN_SUBJECT_COVER` (default `0.0`, off) tests osmotool's
`--min_subject_cover` filter, which discards hits that only cover a small
fraction of the matched reference protein — the failure mode found while
benchmarking proX, where reads matched only the fused-in domain of a much
longer multi-domain protein. Off by default because it's a blunt
instrument (a short read is always much shorter than most full-length
references), so start low (`0.1`-`0.3`) and compare against
`results/metrics/summary.tsv` from a `DIAMOND_MIN_SUBJECT_COVER=0.0` run
to see the recall cost before raising it further.

**Note**: this release also changes osmotool's DIAMOND output format (adds
a `scovhsp` column) — `compute_metrics.py` now expects that column, so
`*.blastx.tsv` files from an older release won't parse correctly against
the updated script. Rerun `benchmark` (not just `all`) against any release
built before this change to regenerate them.

## Adding a new gene family

Capped fetches (every family's negatives, and positives with
`max_positive_override`) are sampled from the whole matching UniProt
population. The accession lists are cached per UniProt release in
`cache/uniprot/`, so the first build after a new UniProt release spends
about an hour listing them and later builds reuse the cache. Pass
`--no-cache` to `fetch_refs.py` to force a fresh listing.

Edit `families.yaml` — add one entry with a `name`, `positive_query`
(UniProt REST query) and `negative_query` (a hard-negative Pfam-family
query). Then rerun `bash run_pipeline.sh all <release_name>`. No other file
needs to change.

See `families.yaml` for query syntax notes (in particular: UniProt's REST
API does **not** support free-text `family:"..."` queries — use
`xref:pfam-PFxxxxx` instead, verified via the
[InterPro/Pfam API](https://www.ebi.ac.uk/interpro/api/entry/pfam/)).

Optional `scope: annotate_only` field: for a family that should still be
built and searchable (e.g. for genome-level co-occurrence checks with
`osmotool annotate`) but excluded from `osmotool profile`'s reported
output — e.g. murB, a near-universal housekeeping gene with no
osmoadaptation-specific signal at the read level. Omit the field (default)
for a normal family available in both modes.

Optional `decoy_from_negatives: true` field: for a family whose
`negative_query` targets one or two *specific named* confusable paralogs
(not a broad Pfam pool) that a single score threshold can't reliably
separate from the true positives — see betL vs betT/caiT, where calibration
proved the two score distributions overlap almost completely. Adds that
family's held-out negative sequences to the DIAMOND db as
`<family>_decoy` references, purely so they can win DIAMOND's existing
best-hit contest away from a mislabeled call; never appears in reported
output. See `modules/local/build_diamond/templates/build_decoy_refs.py`'s docstring for the full
mechanism.

Optional `trim_gt` field (float, default 0.8): per-family override for
trimAl's `-gt` gap threshold in alignment trimming. Worth trying for a
family with a broad, promiscuous Pfam domain (see ectA) where the default
trim may cut away the columns that actually carry the discriminating
signal — validate with the subset-testing workflow below before keeping
it, since it doesn't always help (see ectB's `families.yaml` entry, where
the same experiment made DIAMOND meaningfully worse and was reverted).

> **Status: no family currently uses `pfam_model`.** The pipeline still
> supports it, but all six adoptions (mrpB, mrpE, mrpF, mrpG, gshB, cspA)
> were reverted to custom-built, locally F1-calibrated HMMs on
> 2026-08-01 to 2026-08-04. Pfam's GA cutoffs are tuned on full-length
> proteins and lost recall on short simulated read fragments (gshB HMM F1
> 0.466→0.880 after removal). The trkH/ktrB/ktrD case described below was
> reverted earlier, on 2026-07-24. The guidance below still applies if
> you adopt it for a future family.

Optional `pfam_model: PFxxxxx` field: use that Pfam family's own curated
HMM directly — built from Pfam's hand-curated seed alignment, carrying its
own GA/TC/NC gathering-threshold lines — instead of `hmmbuild`-ing a fresh
model from this repo's own fetched-and-aligned UniProt sequences. Only
worth setting when a gene's real biological target is essentially
co-extensive with a single, specific Pfam family — check via the
[InterPro API](https://www.ebi.ac.uk/interpro/api/entry/pfam/PFxxxxx): a
low `domain_architectures` count and a family name/description that
matches the gene one-to-one are good signs; a domain shared by many
unrelated protein families (e.g. a generic ABC-transporter ATPase or
permease fold) is not, however narrow it looks from protein *count* alone.
Also check no *other* family already in this panel is built from the same
Pfam accession as a genuinely distinct true positive — see `families.yaml`'s
opuAC/opuBC/opuCC vs. proX (all four share PF04069 but are real, separate
paralogs; adopting Pfam's one model+cutoff there would make them
indistinguishable from each other, so none of them use `pfam_model` despite
PF04069 being a reasonably narrow domain by domain-architecture-count
standards). **This isn't just a specificity nicety — for HMM detection
specifically it's a hard failure mode, confirmed in production**: an
earlier version of this panel gave `pfam_model: PF02386` to all three of
trkH/ktrB/ktrD (real, distinct siblings sharing that Pfam family, same
situation as proX/opuAC/opuBC/opuCC). Adopting the same accession for more
than one of them made their HMMs byte-identical — every read scores
exactly the same against all three, `HMM_CALIBRATE`'s alphabetical
`cat hmms/*.hmm` plus `compute_metrics.py`'s best-hit tie-break means
whichever sorts first absorbs 100% of the shared signal, and the rest get
exactly zero HMM recall. The v6 benchmark caught this directly: ktrB
scored real hits, trkH and ktrD scored 0 true positives each, despite
normal, non-trivial DIAMOND recall on the same reads (DIAMOND discriminates
via each family's own distinct reference sequences, not the shared
profile, so it's unaffected). `fetch_refs.py` now fails fast
(`check_no_duplicate_pfam_models`) if two families.yaml entries share a
`pfam_model` accession, specifically to catch this mistake before a full
build wastes time on it again. `PFAM_FETCH` fetches the actual
`.hmm` file (gzipped) from
`https://www.ebi.ac.uk/interpro/api/entry/pfam/<PFxxxxx>/?annotation=hmm` —
confirmed working on that host; note this is **not** the same thing as the
InterPro entry metadata endpoint's `entry_annotations.hmm` counter, which
reports `0` regardless of whether the model is actually fetchable, so don't
use that counter to decide. DIAMOND references are unaffected either way —
those always come from this repo's own UniProt fetch; `pfam_model` only
changes where the HMM comes from. `calibrate_cutoffs.py` reads (and
leaves untouched) the GA line Pfam shipped, rather than computing one, but
still scores this repo's own held-out positive/negative sets against it and
flags the family (`pfam_ga_review_needed` in `cutoff_manifest.tsv`) if
Pfam's cutoff doesn't cleanly separate them — that's a signal the
gene-specificity assumption above was wrong for this particular family, not
that the number needs recalculating; remove `pfam_model` and let it fall
back to a normal custom-built model instead.

**`domain_architectures` count is a screening heuristic, not sufficient on
its own — `qc_scorecard.tsv` review after the first build is a required
step before trusting a `pfam_model` adoption, not a follow-up nicety.**
Confirmed in production: 5 of the 6 `pfam_model` families once adopted in this
panel (mrpB, mrpE, mrpF, mrpG, gshB) got flagged `pfam_ga_review_needed`
despite each passing the InterPro domain-architecture screen beforehand —
only cspA came back clean (`pfam_ga_clean`). A low `domain_architectures`
count and a name/description match are necessary checks, but this repo's
own held-out benchmark is the only thing that actually confirms Pfam's
cutoff behaves well on *this* project's specific negative pool — a
`pfam_model` adoption isn't done until that scorecard has been read.
`qc_scorecard.py` prints an extra, hard-to-miss banner specifically
for `pfam_model` families flagged this way (separate from the general
per-family review list), precisely so this can't be missed by skimming
past a long scorecard run.

Optional `fusion_partner: <other family name>` + `fusion_marker_pfam:
PFxxxxx` fields: for a pair of families whose gene products occur as a
single fused ORF in some lineages instead of two separate genes (currently
only mrpA/mrpB — the Mrp/Mnh Na+/H+ antiporter's A and B subunits, confirmed
via UniProt domain-architecture evidence across Actinomycetota, Bacillota
specifically Paenibacillaceae, and some Alphaproteobacteria — see
`families.yaml`'s mrpA/mrpB entries for the taxonomic detail). Without this,
`check_length_outliers.py` would drop a genuine fused-ORF sequence as a
probable fusion *artifact* (it's roughly the sum of both subunits' lengths)
— exactly wrong here, since Task 1b explicitly wants the fused form
recognized as a real detection target. With `fusion_partner` declared, an
over-length candidate is checked against a second length window centered on
(this family's median + the partner's median) instead of being discarded
outright; `fusion_marker_pfam` (a Pfam domain unique to one partner, absent
from a standalone copy of the other — confirmed via UniProt's own domain
annotation, fetched by `fetch_refs.py` for `fusion_partner` families
specifically) confirms it's a true fusion rather than an unrelated length
outlier when domain evidence is available, and is advisory (not required)
otherwise. Confirmed fusion candidates are routed to
`refs/<family>.fusion_candidates.faa` — excluded from this family's own
MAFFT alignment/HMM (a ~1000aa fused sequence would badly gap the
alignment of ~800aa standalone sequences) but picked up by
`modules/local/build_diamond/templates/build_fusion_refs.py`, which merges a pair's fusion
candidates (deduped by UniProt accession) into one `<familyA>_<familyB>_fused`
DIAMOND-only reference set. The same file also carries fused ORFs pulled
out of a family's own NEGATIVE fetch under a bare locus tag — a fused ORF
matches a family's negative_query just as easily as a true hard negative
does, and can dominate that pool badly enough to distort
`check_length_outliers.py`'s own median-based length filter (confirmed
in production for mrpB: ~68% of its negative pool was genuine mrpA+mrpB
fusion sequences, which became the filter's own "normal" length center and
caused it to discard the true hard negatives instead) — see that script's
docstring for the confirmed case and the domain-evidence-based fix (no
length-based fallback exists for this side, unlike the positive side).
Unlike a decoy reference, a fused-ORF hit is a
**real, reportable** detection target, not a sink for mislabeled calls — it
is not added to `osmotool profile`/`annotate`'s exclusion lists. What a
`_fused` hit should imply for the two individual family calls is
complex-aware scoring logic that lives downstream in `osmotool`, not in
this repo. See `modules/local/release_refs/templates/check_length_outliers.py` and
`modules/local/build_diamond/templates/build_fusion_refs.py` docstrings for the full mechanism.

Optional `paralog_group: <name>` field: close paralogs in the panel that
are expected to confuse each other at read level (e.g. the opu ATPases
opuAA/opuBA/opuCA) share a group name. It changes nothing in the build;
`compute_metrics.py` writes `summary_by_group.tsv`, scoring each group's
members as one unit, so a result can be read per gene or per system.

`report_level: gene | group | none` field (every family in
`families.yaml` sets it): the decision table for `osmotool profile`'s unit
table. `group` sums a paralog group into one row and must be set on all of
its members; otherwise each member is `gene` or `none` on its own.
`scope: annotate_only` requires `none`. The values were set from the
calibration rule in `compute_metrics.py` (`decide_units`): a unit is
reported at the finest level whose `f1_all` on the calibration reads
reaches 0.5 (gene before group), each level taking the better of cascade
off/on, and the cascade config follows the chosen level. `f1_all` counts
every non-background read called as the unit as a false positive, so gene
and group level share one denominator (per-gene `precision` leaves out a
sibling's hard negatives that group precision counts). Step 11c writes the
rule's level next to this field in `osmo_refdb.profile_report_levels.tsv`
and warns when they disagree; that manifest also carries the calibration
numbers for every option, so a choice can be checked without rerunning the
pipeline. A family without the field takes the rule's level (markers.yaml,
subset files).

Optional `max_positive_override: <int>` field: per-family cap on
`fetch_refs.py`'s positive-set fetch, overriding the global
`--max-positive` CLI flag for just that family. Needed for a family whose
`positive_query` is deliberately anchored on a broad Pfam accession rather
than a gene symbol (see cspA below) and would otherwise try to fetch every
UniProt member of that Pfam family.

### A family whose positive set is a whole Pfam family, not a gene symbol

Every family above is anchored on a gene symbol (`gene:xxx`). cspA
(`families.yaml`, 2026-07-23) is deliberately different: its
`positive_query` is `xref:pfam-PF00313 AND taxonomy_id:2` — the Pfam
cold-shock-domain family itself — because the goal is to capture every csp
paralog in a genome (cspA through cspI-type genes and organism-specific
equivalents), not just literally-"cspA"-named orthologs. Pair this with
`max_positive_override` (PF00313 has on the order of 80,000 UniProt member
proteins) if the same pattern is useful for a future family. cspA also
used `pfam_model: PF00313` (Pfam's own CSD gathering cutoff) until
2026-08-04. It was removed because the cutoff cost HMM recall on short
reads (HMM recall 0.627 vs. 0.681 for DIAMOND), so don't copy that part.

### Structural negative-pool contamination: when more data can't fix it

For several families — mrpB, mrpC, mrpE, gshA, gshB, gshF, otsA, mazG,
mscL, and (confirmed later, see below) trkH, ktrB, ktrD — a large
fraction of the "negative" pool drawn from the shared Pfam accession
turns out to be real orthologs of the target gene that UniProt
never assigned a curated gene symbol to (bare locus tags, confirmed by hand
for each family via `refs/<family>.negative.flagged.faa`). This is a
different failure mode from a missing gene-symbol alias (see mrpF/mrpG's
PhaF/yufB fix) or a missing numbered-paralog variant (see murB1/murB2 etc.,
issue #8): those are fixable by widening the query; this one is not,
because there's no alias or symbol to add — the sequences are simply
unannotated at the gene-symbol level.

Confirmed not fixable by adding more data, in either direction:

- Negative side: raising mrpF/mrpG's negative-fetch cap 1000→4000
  (`max_negative_override`) gave 4-5x more *clean* negatives after purity
  filtering but made benchmark F1 measurably worse (mrpF DIAMOND F1
  0.329→0.109), not better — more volume from the same narrow Pfam pool
  dilutes with borderline/confusable cases rather than adding real
  discriminative signal.
- Positive side: the same principle showed up fixing trkA's
  numbered-paralog gap (issue #8) — the added trkA1/trkA2 sequences are
  real orthologs, not an artifact, but their taxonomic diversity widened
  the DIAMOND reference cloud and cost specificity (F1 0.920→0.879).

*Superseded (2026-08-01 to 2026-08-04): the reasoning in this paragraph
did not hold up. mrpB and mrpE were fixed by redesigning the negative pool
(random bacterial background) and gshB by dropping `pfam_model`. All three
now use custom, locally F1-calibrated HMMs. It is kept here as history.*
For `pfam_model` families in this situation (mrpB, mrpE, gshB), this
structural finding *validates* the `pfam_model` choice for a second,
independent reason beyond the domain-architecture screening that
originally justified it: since a locally-calibrated negative set can't be
trusted here, deferring to Pfam's own GA cutoff (trained on Pfam's own
curated seed alignment, not this repo's contaminated pool) is the reliable
path, not just a shortcut. Confirmed via `qc_scorecard.tsv`: all three
carry `hmm_status=pfam_ga_review_needed` rather than
`overlapping_distributions_f1_calibrated`, meaning their HMM path already
sidesteps this problem the same way. (`pfam_ga_review_needed` itself just
flags that Pfam's cutoff needs a human check against this repo's own
held-out sets — see the `pfam_model` section above — not that the model is
broken.)

For families without a `pfam_model` fallback (otsA, mrpC, mazG, mscL,
gshA, and — confirmed via a later re-check, see below — trkH, ktrB,
ktrD), `osmo_refdb.profile_cascade.tsv`'s DIAMOND→HMM fallback does
**not** rescue them: all carry `hmm_status=overlapping_distributions_f1_calibrated`,
meaning the HMM cutoff was calibrated against the exact same contaminated
negative pool as DIAMOND, not an independent clean source. Checked
individually rather than assumed uniformly unreliable, though:

- otsA (DIAMOND precision 0.543) and mrpC (0.479, found during this check
  to be worse than otsA despite not being in the original four-family
  list) — genuinely too low for confident profile-mode reporting.
  `scope: annotate_only` set (same treatment as murB): still built and
  searchable for `osmotool annotate` co-occurrence checks, excluded from
  `osmotool profile`'s reported output.
- mazG: DIAMOND precision 0.870, borderline — kept in profile mode;
  already carries its own caveat in `families.yaml` ("treat any mazG
  signal as suggestive... rather than confirmed osmotic-stress-specific
  function").
- mscL, gshA: DIAMOND precision 0.927 / 0.905 — fine in practice despite
  the contaminated pool; calibration still found a workable separating
  cutoff. No scope change.
- **trkH, ktrB, ktrD** (added after re-checking whether the numbered-paralog
  fix or a decoy conversion could help this group, per `docs/CHANGELOG.md`'s
  Phase 5): same structural contamination pattern, and neither of the two
  usual fixes applies here — the numbered-paralog additions barely moved
  their purity-flag rates (e.g. trkH 27%→29%), and `decoy_from_negatives`
  is specifically ruled out for this trio, since they share PF02386 the
  same close-paralog way proX/opuAC/opuBC/opuCC shared PF04069, and that
  exact clique shape was shown to catastrophically backfire (see the decoy
  section above). No `pfam_model` fallback is possible either (that's
  what removing it fixed the byte-identical-HMM bug for, above). Both
  DIAMOND (0.265–0.293) and HMM (0.121–0.288) precision are low across
  all three, with neither method offering a reliable alternative to the
  other. `scope: annotate_only` set on all three, same criteria as
  otsA/mrpC.

### Test-driving a new family before the full rebuild

A full build and benchmark across every family takes hours. Before committing to that, test-drive
just the families you added: write a `families.yaml` with only those families (and the `hmm_decoys`
they need), and run the workflow on it under a throwaway release name:

```bash
nextflow run main.nf --families_file families_test.yaml --release v3-test
```

It runs the whole build (fetch, QC, calibrate, benchmark) for those families only, and publishes to its
own `releases/v3-test/`. Once you are happy with the result, run the full `families.yaml` under the real
release name.

## Pipeline steps (`pipeline/`)

| Step | Script | Purpose |
|---|---|---|
| 1 | `fetch_refs.py` | Fetch positive + hard-negative sequences per family from UniProt |
| 1b | `check_negative_purity.py` | Flag/drop hard negatives that are actually true positives under a different gene symbol (DIAMOND identity check against the positive set; see `refs/negative_purity_manifest.tsv`) |
| 1c | `check_length_outliers.py` | Flag/drop sequences (positive or negative) whose length is way off their family's median -- likely a multi-domain fusion protein or a partial fragment, not a clean single-domain family member (see `refs/length_outlier_manifest.tsv`) |
| 1d | *retired* | Merged local study sequences (`extra_sequences/`, gitignored, so not reproducible). Removed 2026-11-17 |
| 1e | `add_refseq_positives.py` | For families opted in via `families.yaml`'s `refseq_gene_symbols` (only where UniProt's own gene-symbol population is already fully exhausted and a per-family novelty/contamination check confirmed it's safe), merge additional real positives from NCBI RefSeq into the fetched positive set (see `refs/refseq_positives_manifest.tsv`) |
| 1f | `build_random_background.py` | Sample a shared random bacterial background (default 5000 proteins carrying none of the panel's Pfam domains, minus any >= 70% identical to a panel positive) into `refs/random_background.faa`. Benchmark-only: never used for calibration or in either db. Step 9 simulates reads from it and step 11 reports each family's off-target false positives per million of those reads (see `refs/random_background_manifest.tsv`) |
| 2 | `REFS_CLUSTER` (workflow) | CD-HIT cluster each family's positives and negatives at 90% identity, one set per task (remove redundancy, avoid train/test leakage from near-duplicates); sets of fewer than 5 sequences are kept as they are |
| 3 | `split_train_test.py` | Split positives three ways: train (builds the HMMs and DIAMOND db), calib (sets the HMM GA cutoffs, the DIAMOND bitscore cutoffs and, as simulated reads, the profile cascade config) and test (benchmark only). The held-out 20% (`--test-fraction`) is drawn as before, so train is unchanged from a two-way split; half of it (`--calib-fraction`) becomes calib. Up to v9 there was no calib set and the cutoffs were tuned on the benchmark's own positives. Also splits every family's negatives into train (cutoff calibration in steps 5-6 and 8b; decoy refs in step 8a) / test (benchmark negative reads), so the benchmark never scores negatives that set a cutoff. Non-decoy families hold out half their negatives (`--negative-test-fraction`, default 0.5); decoy families use `--test-fraction` |
| 4 | `MAFFT_ALIGN`, `TRIMAL` (workflow) | MAFFT align + trimAl trim the TRAIN positives, one family per task (fewer than 2 sequences: no alignment). trimAl's gap threshold defaults to 0.8 but can be overridden per family via `families.yaml`'s `trim_gt` |
| 5 | `HMMER_HMMBUILD`, `PFAM_FETCH`, `HMMER_HMMSEARCH` (workflow) | hmmbuild per family (or, for a family marked `pfam_model`, fetch that curated Pfam-A HMM directly instead); score the positive-calib and negative-train sets. HMM decoys: `HMM_DECOYSELECT` (`build_hmm_decoys.py --select-only`), then `CDHIT`, `MAFFT_ALIGN`, `TRIMAL`, `HMMER_HMMBUILD` |
| 6 | `calibrate_cutoffs.py` (`HMM_CALIBRATE`) | Set per-family HMM GA (gathering) cutoffs (or, for `pfam_model` families, read and keep Pfam's own GA line, only flagging it if it doesn't separate this repo's own curated sets well). Where positives and negatives overlap, the cutoff is the median over 1000 seeded bootstrap resamples of the F1-optimal cutoff (`--cutoff-method sweep` for the plain F1 sweep): it moves less between builds. The manifest also records the plain sweep's cutoff, the bootstrap 10th-90th percentile range and `cutoff_range_wide` (the data do not pin the cutoff down) |
| 6b | `qc_scorecard.py` | Consolidate negative-purity, length-outlier, and HMM-cutoff-overlap manifests into one per-family scorecard (`qc_scorecard.tsv`) flagging families that need review |
| 7 | `HMM_CALIBRATE` | Concatenate + hmmpress into one binary HMM database |
| 8a | `build_decoy_refs.py` | For families marked `decoy_from_negatives: true` (e.g. betL), turn their held-out negative-train split into searchable `<family>_decoy` DIAMOND references -- confusable paralogs (e.g. betT/caiT) that can win DIAMOND's best-hit contest away from a mislabeled call, since a single score threshold can't separate them (see betL's `families.yaml` entry) |
| 8d | `build_fusion_refs.py` | For families.yaml pairs marked `fusion_partner` (e.g. mrpA/mrpB), merge each pair's fusion-candidate sequences (flagged by 01c, deduped by UniProt accession) into one searchable `<familyA>_<familyB>_fused` DIAMOND reference -- a real, reportable detection target for a single-ORF fused lineage, not excluded from output the way decoys are |
| 8 | `DIAMOND_TRAINREFS / DIAMOND_MAKEDB` | Build a DIAMOND db from TRAIN positives + any decoy refs from step 8a + any fused refs from step 8d -- excludes `_refseq`-tagged sequences specifically (kept for HMM's alignment, dropped for DIAMOND): confirmed across 4 families that the merged pool helps HMM but mostly hurts DIAMOND's identity-based matching |
| 8b | `calibrate_diamond_cutoffs.py` | Calibrate a per-family minimum DIAMOND bitscore, mirroring HMM's GA cutoff -- DIAMOND otherwise applies one flat `--min_identity` across every family, which can't separate a true gene from a genuinely close paralog (e.g. betL vs betT) whose identity to it happens to exceed that flat threshold |
| 8c | `DIAMOND_SCOPE` | Write `osmo_refdb.profile_excluded_families.txt` (families marked `scope: annotate_only`, e.g. murB, plus all decoy labels) and `osmo_refdb.annotate_excluded_families.txt` (decoy labels only), consumed by `osmotool profile`/`annotate --exclude_families` |
| 9c-11c | steps 9-11 with `--split calib` | Calibration pass, run first by `benchmark`: simulate reads from positive.calib + negative.train, run DIAMOND with `CASCADE=score` (every family flagged at a threshold nothing fails, so no call is dropped but each one's HMM score is kept) and HMM on them. 11c then decides per family or paralog group, on these reads only, the report level and the cascade together: it replays cascade off and on for every unit, applies the report-level rule (see `report_level` above) and writes `osmo_refdb.profile_cascade.tsv` and `osmo_refdb.profile_report_levels.tsv`. Neither was tuned on benchmark reads (the cascade was, up to v9). Output in `results/calib/` |
| 9 | `simulate_reads.py` | Simulate reads from held-out TEST positives + TEST negatives (the negative split disjoint from calibration and from any decoy refs, so a benchmark read is never scored against a threshold or db entry built from its own source sequence) |
| 9b | `compute_read_truth.py` | Locate each read on its source contig (k-mer match) and label it by actual CDS overlap, not just its source construct -- a read that lands entirely in flanking DNA carries no gene signal and must not be scored as a positive. Also labels real background reads and records read length. |
| 10 | `run_benchmark.sh` | Run `osmotool profile` with its default DIAMOND+HMM cascade (`CASCADE=off` for DIAMOND alone, `CASCADE=score` for the calibration pass; `SKIP_HMM=1` reruns only the DIAMOND side) and HMM (hmmscan; `JOBS` samples in parallel) on the simulated reads. `--keep_aln` keeps the raw DIAMOND alignments and the cascade's hmmsearch table, so one run scores both DIAMOND alone and DIAMOND+cascade |
| 11 | `compute_metrics.py` | Precision/recall/F1, ROC/PR curves, per-family best-threshold summary, read-length-stratified summary and background false-positive rate, for three methods: `diamond`, `diamond_cascade` (what `osmotool profile` reports by default; per-read calls rebuilt from osmotool's kept files and checked against its `gene_counts.tsv`) and `hmm`. In the calibration pass (no `--cascade-config`) it instead decides and writes `osmo_refdb.profile_cascade.tsv` and `osmo_refdb.profile_report_levels.tsv` (step 9c-11c); on the test reads it scores those choices per reporting unit in `summary_by_report_unit.tsv` (the numbers for the paper) and checks osmotool's `unit_counts.tsv` against the rebuilt calls. `summary.tsv` counts false positives three ways: `fp` (the family's own hard negatives and flanks; `precision`/`f1`, comparable to releases up to v8), `fp_cross` (other panel genes' real reads called as this family, mostly close-paralog confusion; `precision_panel`/`f1_panel`), and `fp_background`/`background_fp_per_million` (random background reads). `fp_other_negatives` (other families' hard negatives called as this family) is in none of those, since sibling families share one negative pool. `fp_all` = all three, with `precision_all`/`f1_all`: the same denominator at gene and group level, used for the report-level and cascade decision. `hard_negative_fp_per_million` and `background_fp_per_million` are rates per million negative reads, independent of how many negatives were mixed in: use those to compare families or releases. Families sharing a `families.yaml` `paralog_group` are also scored as one unit in `summary_by_group.tsv`, where a read called as a sibling counts as correct (`within_group_swaps`) |
| 11b | `modules/local/benchmark/templates/benchmark_annotate.py` | Annotate-mode benchmark at protein level: runs `osmotool annotate --proteins` (HMM with GA cutoffs, and DIAMOND with its bitscore cutoffs) on the held-out test positives, test hard negatives and random background, and scores per-protein calls from `gene_coordinates.tsv`. Output in `results/annotate/` (`summary.tsv`, `summary_by_group.tsv`, `overall.tsv`) |

| 12 | `MARKERS` (workflow) | Build the marker genes (`markers.yaml`, see below) apart from the osmo families: their own references, split and DIAMOND database, packaged as the release's `markers/` (`package_markers.py`): DIAMOND db, its source FASTA (for reference lengths), excluded decoy labels, a per-marker manifest |
| 14 | `package_catalog.py` (`run_pipeline.sh catalog`) | Build a gene catalog into the release as `catalogs/<name>/` with `osmotool build-catalog`: translate the catalog genes, annotate them with this release's HMMs under their gathering cutoffs (as `annotate --method hmm`), index the whole catalog with bwa-mem3. Adds provenance (source, input checksums, release) to `catalog.json`, the gzipped hmmsearch domain table, and `SHA256SUMS`. The catalog is part of the release: rebuild it for every release |

### Gene catalogs

(Running this on a cluster with a full-size gut catalog: `docs/hpc.md`.
The prevalence-filtered GMGC catalogs are on Zenodo, record 10473258;
`modules/local/catalog/fetch/templates/fetch_zenodo_member.py` downloads one habitat without the 19.3 GB zip.)

`osmotool profile --catalog NAME` assigns reads to a gene catalog the
release ships under `catalogs/NAME/` before translated search, the approach
of catalog-based profilers such as Cayman: the family of each catalog gene is
decided once, at full length, and reads are then matched at >= 97%
nucleotide identity. The catalog's annotation is built into the release
and versioned with it: it depends on this release's HMMs and cutoffs, so
users never annotate anything themselves.

```bash
CATALOG_NAME=gmgc_gut \
CATALOG_GENES=/data/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz \
CATALOG_SOURCE="GMGC v1.0, human gut, prevalence > 0.5%" \
CATALOG_SOURCE_URL=https://zenodo.org/records/10473258 THREADS=32 \
bash run_pipeline.sh catalog v10
```

- The release ships an index of only the genes with an HMM hit
  (`CATALOG_INDEX=hits`, the default), so it stays small and the catalog
  works straight away. Sub-cutoff homologs and decoys are in it, so their
  reads are recognised and dropped; reads from genes without any hit fall
  back to translated search.
- The full, competitive index of every catalog gene (bwa-mem3) does not
  depend on the release, only on the gene file and the aligner. It is about
  3.4 bytes per byte of FASTA on the catalogs measured so far (a 137 MB
  sample of the human-gut catalog gave a 472 MB index); the size for a whole
  gut catalog has not been measured (see `docs/hpc.md`). It is built separately
  and distributed on its own:

  ```bash
  osmotool index-catalog --genes GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz \
      -o gmgc_gut_index
  osmotool profile releases/v10 -1 R1.fq.gz -2 R2.fq.gz \
      --catalog gmgc_gut --catalog_index gmgc_gut_index -o sample
  ```

  osmotool checks that the index was built from the same gene file the
  release's catalog was annotated from (`genes_sha256`).
- The catalog is not part of the simulated-read benchmark yet, and
  `--catalog` is not osmotool's default.

Steps 4, 5 (hmmbuild) and 8 only touch the TRAIN split. Steps 5-6 and 8b
(cutoffs) and the calibration pass 9c-11c (cascade config) only touch the
CALIB positives and TRAIN negatives. Steps 9-11b only touch the TEST split.
Nothing that is tuned sees a sequence the benchmark scores.

`REUSE_REFS_FROM=releases/<name>` rebuilds from another release's fetched
and clustered sequences (skipping steps 1-2) and, unless
`REUSE_ALIGNMENTS=0`, its alignments. `releases/v10` was built this way from
v9: same sequences, same HMMs and DIAMOND db, new cutoffs and cascade config.
Reusing the alignments matters because multi-threaded MAFFT is not
reproducible from run to run.

## Marker genes for copies-per-genome normalisation

`osmotool profile` reports RPM, which shifts with host/eukaryotic DNA
content and average genome size across samples. When a release ships a
`markers/` subdirectory, `profile` also reports `copies_per_genome`: a
family's reads per kb of effective gene length divided by the mean of ten
universal single-copy marker genes, i.e. the average number of copies per
genome in the community. The method follows MicrobeCensus (Nayfach &
Pollard 2015), reimplemented for DIAMOND rather than vendored (it's GPL-3.0,
and its trained cutoffs/coefficients only apply to RAPsearch2 scores); see
`osmotool`'s `normalise.py` for the details.

The markers (`markers.yaml`: the ten fetchMG/mOTUs COGs -- ychF, pheS,
argS, serS, cysS, leuS, valS, tsaD, ffh, ftsY) go through the same fetch ->
QC -> cluster -> split -> calibrate -> benchmark pipeline as the osmo
families, but as **their own release with their own DIAMOND db**, so they
never compete with osmo families in best-hit selection and adding them
can't change any osmo family's benchmark. murB is deliberately not a
marker: it has a reported salt-tolerance association (Culligan et al.
2012) and real paralogs in some taxa, and a normaliser must not covary
with what it normalises.

The workflow builds them with every release (`MARKERS`): `markers.yaml`'s references, split and DIAMOND
database, apart from the osmo families, packaged as `releases/<release>/markers/` and given to every
osmotool step after the gate. The markers have no benchmark of their own in the workflow, so the marker
manifest's DIAMOND precision/recall columns are empty: a marker with poor DIAMOND precision would inflate
the denominator for every family, which the old separate marker build's benchmark used to check.

Phase-1 limitation: per-family detection sensitivity isn't corrected for.
Markers and osmo families use the same DIAMOND thresholds, so the bias
partly cancels, but `copies_per_genome` is best read as comparable across
samples rather than as an exact absolute copy number. Calibrating
per-family sensitivity coefficients from the simulated-read benchmark
(MicrobeCensus's approach) is the planned phase 2.

## Real-genome validation

The `benchmark` step above measures performance on reads simulated from
UniProt-derived sequences. `docs/real_genome_validation.md` holds a
separate, complementary (semi-manual) workflow for checking gene calls
against a real, independently-annotated genome instead — see
`docs/real_genome_validation.md` for the full procedure.

## Releasing a database for `osmotool` to use

Once you're happy with a release's benchmark results, publish
`releases/<name>/osmo_refdb.dmnd` and `releases/<name>/hmms/osmo_refdb.hmm`
(+ its `.h3*` index files) as a versioned data release (e.g. a GitHub
release asset or Zenodo record) that `osmotool` users download and point
`--database` at. Cite the specific release version in any publication using
`osmotool`.

**Latest release: v5** — DOI: [10.5281/zenodo.21420253](https://doi.org/10.5281/zenodo.21420253).
Download the full release (databases + cutoffs/exclusion files + full
build/benchmark provenance) from that record. Cite this DOI when using
`osmotool` with this database.
