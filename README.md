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

`osmotool profile` (short reads) never needs HMMER -- HMM's calibrated
cutoffs are scaled for full-length sequences, not read fragments, so
DIAMOND alone is used there. `osmotool annotate` (full-length ORFs from
an assembly) defaults to HMM (`--method hmm`) since real-genome testing
showed it meaningfully more specific than DIAMOND for this reference
database (see `hmms/osmo_refdb.hmm` + `--hmm_db`); pass `--method diamond`
if you want to avoid the HMMER dependency there too.

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
a `scovhsp` column) — `11_compute_metrics.py` now expects that column, so
`*.blastx.tsv` files from an older release won't parse correctly against
the updated script. Rerun `benchmark` (not just `all`) against any release
built before this change to regenerate them.

## Adding a new gene family

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
output. See `pipeline/08a_build_decoy_refs.py`'s docstring for the full
mechanism.

Optional `trim_gt` field (float, default 0.8): per-family override for
trimAl's `-gt` gap threshold in alignment trimming. Worth trying for a
family with a broad, promiscuous Pfam domain (see ectA) where the default
trim may cut away the columns that actually carry the discriminating
signal — validate with the subset-testing workflow below before keeping
it, since it doesn't always help (see ectB's `families.yaml` entry, where
the same experiment made DIAMOND meaningfully worse and was reverted).

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
exactly the same against all three, `07_press_hmms.sh`'s alphabetical
`cat hmms/*.hmm` plus `11_compute_metrics.py`'s best-hit tie-break means
whichever sorts first absorbs 100% of the shared signal, and the rest get
exactly zero HMM recall. The v6 benchmark caught this directly: ktrB
scored real hits, trkH and ktrD scored 0 true positives each, despite
normal, non-trivial DIAMOND recall on the same reads (DIAMOND discriminates
via each family's own distinct reference sequences, not the shared
profile, so it's unaffected). `01_fetch_refs.py` now fails fast
(`check_no_duplicate_pfam_models`) if two families.yaml entries share a
`pfam_model` accession, specifically to catch this mistake before a full
build wastes time on it again. `05_build_hmms.sh` fetches the actual
`.hmm` file (gzipped) from
`https://www.ebi.ac.uk/interpro/api/entry/pfam/<PFxxxxx>/?annotation=hmm` —
confirmed working on that host; note this is **not** the same thing as the
InterPro entry metadata endpoint's `entry_annotations.hmm` counter, which
reports `0` regardless of whether the model is actually fetchable, so don't
use that counter to decide. DIAMOND references are unaffected either way —
those always come from this repo's own UniProt fetch; `pfam_model` only
changes where the HMM comes from. `06_calibrate_cutoffs.py` reads (and
leaves untouched) the GA line Pfam shipped, rather than computing one, but
still scores this repo's own held-out positive/negative sets against it and
flags the family (`pfam_ga_review_needed` in `cutoff_manifest.tsv`) if
Pfam's cutoff doesn't cleanly separate them — that's a signal the
gene-specificity assumption above was wrong for this particular family, not
that the number needs recalculating; remove `pfam_model` and let it fall
back to a normal custom-built model instead.

Optional `fusion_partner: <other family name>` + `fusion_marker_pfam:
PFxxxxx` fields: for a pair of families whose gene products occur as a
single fused ORF in some lineages instead of two separate genes (currently
only mrpA/mrpD — Bacillota's Mrp/Mnh Na+/H+ antiporter subunits A and D).
Without this, `01c_check_length_outliers.py` would drop a genuine fused-ORF
sequence as a probable fusion *artifact* (it's roughly the sum of both
subunits' lengths) — exactly wrong here, since Task 1b explicitly wants the
fused form recognized as a real detection target. With `fusion_partner`
declared, an over-length candidate is checked against a second length
window centered on (this family's median + the partner's median) instead of
being discarded outright; `fusion_marker_pfam` (a Pfam domain unique to one
partner, absent from a standalone copy of the other — confirmed via
UniProt's own domain annotation, fetched by `01_fetch_refs.py` for
`fusion_partner` families specifically) confirms it's a true fusion rather
than an unrelated length outlier when domain evidence is available, and is
advisory (not required) otherwise. Confirmed fusion candidates are routed
to `refs/<family>.positive.fusion_candidates.faa` — excluded from this
family's own MAFFT alignment/HMM (a ~1300aa fused sequence would badly gap
the alignment of ~800aa standalone sequences) but picked up by
`pipeline/08d_build_fusion_refs.py`, which merges a pair's fusion
candidates (deduped by UniProt accession) into one `<familyA>_<familyB>_fused`
DIAMOND-only reference set. Unlike a decoy reference, a fused-ORF hit is a
**real, reportable** detection target, not a sink for mislabeled calls — it
is not added to `osmotool profile`/`annotate`'s exclusion lists. What a
`_fused` hit should imply for the two individual family calls is
complex-aware scoring logic that lives downstream in `osmotool`, not in
this repo. See `pipeline/01c_check_length_outliers.py` and
`pipeline/08d_build_fusion_refs.py` docstrings for the full mechanism, and
`families.yaml`'s mrpA entry for an honest note on searching UniProt
directly for confirmed mrpA/mrpD fusion examples (none found at time of
writing — the mechanism is built defensively per the spec regardless).

Optional `max_positive_override: <int>` field: per-family cap on
`01_fetch_refs.py`'s positive-set fetch, overriding the global
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
proteins) and `pfam_model` (Pfam's own long-established CSD gathering
cutoff) if the same pattern is useful for a future family.

### Test-driving a new family before the full rebuild

A full `build` + `benchmark` across every family takes on the order of 90
minutes. Before committing to that, test-drive just the family/families
you added:

```bash
python pipeline/make_family_subset.py --families families.yaml \
    --only galE,mazG,murB --out families_test.yaml

FAMILIES_FILE=families_test.yaml bash run_pipeline.sh all v3-test
```

This runs the entire pipeline (fetch, QC, calibrate, benchmark) but only
for the named families, finishing in minutes instead of hours. Use a
throwaway `release_name` (like `v3-test` above, not `v3`) so it writes to
its own `releases/` directory and can never leave your real release in a
partial state. Once you're happy with the result, rerun the full pipeline
without `FAMILIES_FILE` against your real release name to rebuild
everything together.

## Pipeline steps (`pipeline/`)

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_fetch_refs.py` | Fetch positive + hard-negative sequences per family from UniProt |
| 1b | `01b_check_negative_purity.py` | Flag/drop hard negatives that are actually true positives under a different gene symbol (DIAMOND identity check against the positive set; see `refs/negative_purity_manifest.tsv`) |
| 1c | `01c_check_length_outliers.py` | Flag/drop sequences (positive or negative) whose length is way off their family's median -- likely a multi-domain fusion protein or a partial fragment, not a clean single-domain family member (see `refs/length_outlier_manifest.tsv`) |
| 2 | `02_cluster_cdhit.sh` | CD-HIT cluster positives at 90% identity (remove redundancy) |
| 3 | `03_split_train_test.py` | Split positives into train (build) / test (held-out benchmark). For families marked `decoy_from_negatives: true`, also splits their negative set the same way (see step 8a) |
| 4 | `04_align_trim.sh` | MAFFT align + trimAl trim the TRAIN positives. trimAl's gap threshold defaults to 0.8 but can be overridden per family via `families.yaml`'s `trim_gt` |
| 5 | `05_build_hmms.sh` | hmmbuild per family (or, for a family marked `pfam_model`, fetch that curated Pfam-A HMM directly instead); score positive-test/negative sets |
| 6 | `06_calibrate_cutoffs.py` | Set per-family HMM GA (gathering) cutoffs (or, for `pfam_model` families, read and keep Pfam's own GA line, only flagging it if it doesn't separate this repo's own curated sets well) |
| 6b | `06b_qc_scorecard.py` | Consolidate negative-purity, length-outlier, and HMM-cutoff-overlap manifests into one per-family scorecard (`qc_scorecard.tsv`) flagging families that need review |
| 7 | `07_press_hmms.sh` | Concatenate + hmmpress into one binary HMM database |
| 8a | `08a_build_decoy_refs.py` | For families marked `decoy_from_negatives: true` (e.g. betL), turn their held-out negative-train split into searchable `<family>_decoy` DIAMOND references -- confusable paralogs (e.g. betT/caiT) that can win DIAMOND's best-hit contest away from a mislabeled call, since a single score threshold can't separate them (see betL's `families.yaml` entry) |
| 8d | `08d_build_fusion_refs.py` | For families.yaml pairs marked `fusion_partner` (e.g. mrpA/mrpD), merge each pair's fusion-candidate sequences (flagged by 01c, deduped by UniProt accession) into one searchable `<familyA>_<familyB>_fused` DIAMOND reference -- a real, reportable detection target for a single-ORF fused lineage, not excluded from output the way decoys are |
| 8 | `08_build_diamond_db.sh` | Build a DIAMOND db from TRAIN positives + any decoy refs from step 8a + any fused refs from step 8d |
| 8b | `08b_calibrate_diamond_cutoffs.py` | Calibrate a per-family minimum DIAMOND bitscore, mirroring HMM's GA cutoff -- DIAMOND otherwise applies one flat `--min_identity` across every family, which can't separate a true gene from a genuinely close paralog (e.g. betL vs betT) whose identity to it happens to exceed that flat threshold |
| 8c | `08c_write_scope_manifest.py` | Write `osmo_refdb.profile_excluded_families.txt` (families marked `scope: annotate_only`, e.g. murB, plus all decoy labels) and `osmo_refdb.annotate_excluded_families.txt` (decoy labels only), consumed by `osmotool profile`/`annotate --exclude_families` |
| 9 | `09_simulate_reads.py` | Simulate reads from held-out TEST positives + negatives (for decoy families, from the negative-test split disjoint from what became decoy refs, so a benchmark read can never trivially match its own literal source sequence sitting in the db) |
| 9b | `09b_compute_read_truth.py` | Locate each read on its source contig (k-mer match) and label it by actual CDS overlap, not just its source construct -- a read that lands entirely in flanking DNA carries no gene signal and must not be scored as a positive. Also labels real background reads and records read length. |
| 10 | `10_run_benchmark.sh` | Run DIAMOND (`osmotool profile`) + HMM (hmmscan) on simulated reads |
| 11 | `11_compute_metrics.py` | Precision/recall/F1, ROC/PR curves, per-family best-threshold summary, read-length-stratified summary, background false-positive rate, and `osmo_refdb.profile_cascade.tsv` (families where DIAMOND's benchmark precision falls below a threshold, paired with HMM's best short-read bitscore cutoff -- the config `osmotool profile`'s DIAMOND+HMM cascade consumes) |

Steps 4–8 only ever touch the TRAIN split; steps 9–11 only ever touch the
TEST split — this train/test separation is required so the benchmark
measures genuine generalization rather than memorization.

## Real-genome validation

The `benchmark` step above measures performance on reads simulated from
UniProt-derived sequences. `pipeline/real_genome_validation/` holds a
separate, complementary (semi-manual) workflow for checking gene calls
against a real, independently-annotated genome instead — see
`pipeline/real_genome_validation/README.md` for the full procedure.

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
