# Changelog

Chronological summary of osmo_refdb's development, from the initial
pipeline build through the current negative-pool quality follow-up work
on branch `cluster-negatives-cdhit`. Full rationale and evidence for
each item lives in `families.yaml`'s per-family `description` fields and
`README.md`; this file is a summary, not the source of truth.

## v5 release — QC infrastructure and pipeline foundations

Consolidated the project's accumulated pipeline work into its first real
checkpoint (`6f60d1a`):

- **QC infrastructure**: `01b`/`01c` flag/drop hard negatives that are
  secretly true positives (identity check against positives) and length
  outliers (fusion proteins/fragments); `06_calibrate_cutoffs.py` moved
  from a 99th-percentile-of-negatives heuristic to an F1-maximizing
  threshold sweep for HMM GA cutoffs, with a minimum-negative-sample-size
  guard; `06b_qc_scorecard.py` added a consolidated per-family review
  scorecard; `08b_calibrate_diamond_cutoffs.py` added per-family DIAMOND
  bitscore cutoffs (mirroring the HMM side, since a single flat
  `--min_identity` can't separate close paralogs); `09b_compute_read_truth.py`
  moved to CDS-overlap-based read truth instead of blanket per-construct
  labeling; `11_compute_metrics.py` added read-length-stratified metrics
  and the DIAMOND+HMM profile-mode cascade config.
- **New families**: galE, mazG, murB (Culligan et al. 2012 ISME J
  salt-tolerance locus), murB scoped `annotate_only`.
- **betL decoy fix**: betL's positive/negative score distributions
  overlapped almost completely with its named confusable paralogs
  betT/caiT — no threshold could separate them. Introduced
  `decoy_from_negatives`, adding betL's own QC'd negative pool to the
  DIAMOND db as excluded-from-output decoy references so the existing
  best-hit contest could win real betT/caiT reads away from mislabeled
  betL calls. Required splitting negatives into train/test
  (`03_split_train_test.py`) and independent per-family RNG seeding.
  Validated: betL DIAMOND F1 0.291→0.897.
- **ectA/ectB trim_gt experiment**: per-family trimAl gap-threshold
  override. Kept for ectA (F1 0.773→0.806); reverted for ectB after it
  made DIAMOND meaningfully worse — an early example of "tried, measured,
  reverted" discipline this project uses throughout.
- **Tooling**: `make_family_subset.py` + `FAMILIES_FILE` for testing a
  single family in minutes instead of a full ~90-minute rebuild;
  `real_genome_validation/` formalized checking gene calls against real,
  independently-annotated genomes rather than only synthetic
  UniProt-derived benchmarks.

## PR #6 — Task 1 gene panel (Firmicutes gap-filling)

Closed two confirmed gaps: Bacillota lacking a phylum-appropriate
specificity gene, and nhaA being functionally absent in Bacillota.

- **Task 1a** (`1780b9b`) — Firmicutes compatible-solute transporters:
  the *B. subtilis*-type glycine betaine/carnitine/choline ABC
  transporter systems (opuA/opuB/opuC, 9 subunit families). Updated
  proX's `negative_query` to cross-exclude the new opuAC/opuBC/opuCC
  symbols, since all four families share Pfam domain PF04069 as real,
  distinct true positives.
- **Task 1b** (`bbd6972`) — Mrp/Mnh Na+/H+ antiporter complex (mrpA-mrpG):
  Bacillota's functional substitute for nhaA. Introduced `pfam_model`
  (use a Pfam family's own curated HMM+GA cutoff directly, after
  per-gene InterPro specificity checks) and `fusion_partner`/
  `fusion_marker_pfam` (for gene pairs that occur as a single fused ORF
  in some lineages).
- **Task 1c** (`6b27254`) — gut-validated osmotic-stress subsystems:
  betA/betB (choline→glycine betaine oxidation), gshA/gshB/gor
  (glutathione biosynthesis), cspA-family cold shock proteins. Basis: Ng
  et al. 2023 (mBio), in vivo humanized-mouse data. cspA is the one
  family in this panel anchored on the bare Pfam accession
  (`xref:pfam-PF00313`) rather than a gene symbol, to capture every csp
  paralog per genome; introduced `max_positive_override` to cap its
  ~80,000-member fetch.
- **Task 1d** (`7eadb01`) — Trk/Ktr constitutive K+ uptake, an
  alternative route to the inducible Kdp system already covered by
  kdpA. Found a naming trap during research: the literal gene symbol
  "trkA" mostly does NOT refer to this system in Firmicutes — the real
  Firmicutes-lineage ortholog is consistently named ktrA/ktrB/ktrD
  instead.

**Bugs found and fixed while building out Task 1** (all same-day
follow-ups, issues #1-#5):

- `d175912` — **hmmpress SSI collision**: two families sharing a
  `pfam_model` accession (trkH/ktrB/ktrD → PF02386) produced
  byte-identical ACC lines and hmmpress failed outright. Fixed by
  rewriting ACC to the family name, same as NAME.
- `de60d2e` — **trkA/trkH vs. ktrA/ktrB/ktrD split**: originally folded
  ktrA/ktrB/ktrD in as gene-symbol synonyms of trkA/trkH (the betL
  pattern). Corrected to five separate families, matching how nhaA vs.
  mrpA-mrpG are already kept separate for the same antiporter function
  rather than merged.
- `66b100b` — **trkH/ktrB/ktrD zero HMM recall** (issue #2): all three
  adopted the same `pfam_model` accession, making their HMMs
  byte-identical — hmmscan couldn't discriminate between them at all,
  and the alphabetically-first family (ktrB) absorbed 100% of the shared
  signal via the best-hit tie-break. Confirmed in the v6 benchmark: trkH
  and ktrD scored exactly 0 true positives each despite normal DIAMOND
  recall on the same reads. Fixed by reverting to independently
  custom-built HMMs for all three; added `check_no_duplicate_pfam_models()`
  to fail fast if this happens again.
- `5602c06` — **fused-ORF pairing correction** (issue #3): Task 1b's
  spec text named mrpA+mrpD as a fused-ORF pair, from general
  recollection rather than verification, and was wrong. Direct UniProt
  search found no confirmed mrpA-mrpD fusion, but did find a different,
  real one: mrpA+mrpB, dispersed across Actinomycetota, Bacillota
  (specifically Paenibacillaceae), and some Alphaproteobacteria — with
  genus *Bacillus* itself never showing the fused form.
- `2e19d24` — **UniProt fetch order bias** (issue #4): capped fetches
  paginated UniProt in its default result order, which is not
  representative (confirmed: a plain cspA-style query put E. coli/B.
  subtilis in 8 of its first 10 hits) — systematically skewing every
  capped family's reference/calibration set toward heavily-sequenced lab
  organisms. Fixed via bounded oversampling on `sort=accession asc`
  order rather than exhaustive enumeration (which doesn't scale — ectA's
  negative pool alone has 311,116 UniProt members). Added
  `diversity_stats()` reporting so a regression like this shows up in
  the build's own output next time.
- `faeddd6` — **gshF addition** (issue #5): several common gut genera
  (Streptococcus, Enterococcus, Listeria, Clostridium) use a single
  bifunctional gshF/gshAB enzyme instead of separate gshA+gshB genes,
  which without a dedicated family showed a false "no glutathione
  biosynthesis" signal in exactly the genera gshA/gshB's evidence-tier
  claim depends on. Confirmed gshF is a genuinely distinct third enzyme
  family (different synthetase domain than gshB), not a literal fusion
  of gshA+gshB — modeled as an ordinary family, not `fusion_partner`.
- `eb60e97` — **heredoc backtick fix**: a markdown-style backtick in a
  Python comment embedded in a bash heredoc was interpreted as shell
  command substitution (harmless in practice, but a stray error in every
  `pfam_model` family's build log).

## PR #7 — mrpF/mrpG PhaF/yufB synonym fix

`8f1fbd5`: both families' descriptions documented a third gene-symbol
alias from the start (PhaF for mrpF, yufB for mrpG) but never actually
carried it into `positive_query`/`negative_query` — a straightforward
query bug, independent of the `pfam_model` adoption itself.

- mrpF: confirmed real and consequential — 142 PhaF-tagged genuine
  PF04066 entries were sitting in the "negative" pool. Fixing it gave a
  real, confirmed benchmark improvement (DIAMOND F1 0.244→0.329, HMM F1
  0.193→0.250).
- mrpG: added for consistency, but verified this does **not** explain
  mrpG's own 59% purity contamination the way PhaF explains mrpF's —
  every yufB-tagged entry already also carried the mrpG gene name on the
  same record. F1 barely moved, within noise.
- Both families' remaining negative-pool contamination was traced to the
  same root structural cause later formalized in this branch's Phase 3
  (see below): PF04066/PF03334 have few members genuinely distinct from
  mrpF/mrpG themselves.

## This branch (`fix-numbered-paralog-gaps`) — negative-pool quality follow-ups

PR #6's gene panel and PR #7's mrpF/mrpG fix surfaced a family of related
data-quality issues in how per-family hard-negative pools are built. The
sections below cover the investigation and fix work that followed.

### Investigation (Context)

Before any fix landed, evidence was pulled from the last clean `v6` build
(`releases/v6/refs/negative_purity_manifest.tsv`, flagged-sequence gene
lookups via UniProt) and one experiment was run live:

- **`max_negative_override` disproven.** Raising mrpF/mrpG's negative-fetch
  cap 1000→4000 gave 4-5x more *clean* negatives after purity filtering but
  made benchmark F1 measurably **worse** (mrpF DIAMOND F1 0.329→0.109). More
  volume from the same narrow Pfam pool dilutes with borderline/confusable
  cases rather than adding discriminative signal. Change reverted, finding
  documented instead (see Phase 3).
- Flagged-negative sequences for mrpB/mrpC/mrpE, gshB, otsA, mazG, mscL all
  traced to the same root cause: real orthologs UniProt never assigned a
  curated gene symbol to (bare locus tags), not a synonym gap or threshold
  problem — see Phase 3.
- Two real, cheap-to-fix gaps surfaced along the way: **numbered-paralog**
  gene symbols (`murB1`/`murB2`, `trkH1`/`trkH2`/etc.) missing from
  `positive_query` entirely — see Phase 2.
- proX/opuAC/opuBC/opuCC showed low purity-flag rates but poor precision —
  genuine cross-paralog confusion between four real, distinct genes — see
  Phase 4.

### Phase 1 — Cleanup

Discarded the `max_negative_override` experiment's code changes (kept the
finding, not the code). No behavior change.

### Phase 2 — Numbered-paralog gap fix (issue #8, PR #9)

- Started as a 6-family fix (murB, trkA, trkH, ktrA, ktrB, ktrD), each
  numbered variant verified per-accession against UniProt to confirm it
  maps to the correct sibling family.
- Subset-tested against `v6`: **not** a uniform win. murB/trkH improved;
  ktrA/ktrB/ktrD's swings were consistent with noise from small held-out
  test sets; **trkA regressed** (DIAMOND F1 0.920→0.879).
- trkA's regression was root-caused, not just observed: of 372 UniProt
  `trkA1`/`trkA2` entries, 82% (median 223aa) are short single-domain
  fragments already excluded by the existing length-outlier filter — the
  68 that survive are real, taxonomically diverse full-length orthologs
  whose diversity widens the DIAMOND reference cloud and costs specificity.
  Kept anyway (real genes, not an artifact); documented in `families.yaml`.
- **Expanded to a full systematic sweep** of all 43 families after
  discovering this bug class is near-universal, not rare: every
  single-token gene symbol in the panel has real, curated numbered
  variants in UniProt, from a handful (mazG: 9) up to hundreds (nhaA:
  304/389, galE: 678/530). Verified as real curated data via direct
  UniProt lookups, not a search-engine artifact. `cspA` is structurally
  exempt (anchored on the bare Pfam accession, not a gene symbol, by
  design).
- Added `<symbol>1`/`<symbol>2` variants to 36 of 43 families'
  `positive_query`/`negative_query`, including per-accession
  cross-contamination checks for every family sharing a Pfam accession
  with a sibling (proX/opuAC/opuBC/opuCC, opuAA/opuBA/opuCA,
  opuAB/opuBB/opuCB, mrpA/mrpD, gshA/gshF) — all checked clean.
- Given the scope (36/43 families), treated as a full rebuild (`releases/v7`)
  rather than a subset test, per the project convention of using a new
  version number for any full rebuild.

### Phase 3 — Structural negative-pool contamination (documentation)

Documented, for mrpB/mrpC/mrpE/gshA/gshB/gshF/otsA/mazG/mscL, that their
negative-pool contamination is structural (real, unlabeled orthologs from
the same Pfam family) and **not** fixable by more data — same principle as
the Phase 1 `max_negative_override` disproof, and as Phase 2's trkA
finding, just triggered from different directions.

- **`pfam_model` vindication**: for mrpB, mrpE, gshB (all flagged
  `pfam_ga_review_needed` in `cutoff_manifest.tsv`), this structural
  finding validates the original `pfam_model` adoption for a second,
  independent reason — a locally-calibrated negative set couldn't have
  been trusted here anyway, so deferring to Pfam's own GA cutoff is the
  reliable path.
- **Production disposition, checked individually rather than assumed
  uniform**: of the four families with no `pfam_model` fallback (otsA,
  mazG, mscL, gshA), `profile_cascade.tsv`'s DIAMOND→HMM fallback doesn't
  rescue any of them (HMM calibration draws on the same contaminated
  pool). Only **otsA** (DIAMOND precision 0.543) and, on the same
  criteria, **mrpC** (0.479, lower still) needed a scope change — set to
  `scope: annotate_only`. mazG (0.870) and mscL/gshA (0.927/0.905) are
  fine in practice; calibration still found a workable cutoff despite the
  noisy pool.
- Added a new README section, "Structural negative-pool contamination:
  when more data can't fix it," documenting the pattern for future family
  additions.

### Phase 4 — proX/opuAC/opuBC/opuCC decoy conversion (tried, disproven, reverted)

- Redesigned all four families to target each other as decoy sources
  (mirroring betL vs. betT/caiT) instead of drawing negatives from the
  broad, anonymous PF04069 pool, with `decoy_from_negatives: true` set on
  all four. Added a dedup check to `08a_build_decoy_refs.py` first
  (`load_all_positive_accessions()`, kept in the codebase) so a decoy set
  can't contain the exact same sequence as a real positive reference
  elsewhere in the combined DIAMOND db.
- **Subset-tested (`releases/v-test-phase4`) before merging — result was a
  catastrophic DIAMOND recall collapse for all four families**, not the
  modest recall cost anticipated: opuBC recall dropped to exactly 0.0;
  proX/opuAC/opuCC all fell under 0.08. Root-caused via
  `opuBC.positive.gene_counts.tsv`: opuBC's own genuine positive reads
  scored best against the *other* families' decoy entries, not opuBC's
  own thin (93-member) reference set.
- Mechanism: unlike betT/caiT (only superficially score-confusable with
  betL, genuinely distinguishable at the sequence level),
  proX/opuAC/opuBC/opuCC are close enough paralogs of *each other* that
  adding them back into the searchable pool as decoys lets them
  systematically outcompete a family's own thin reference set. The dedup
  fix addressed the specific collision it was designed for (same
  accession as both decoy and positive reference) but not this failure
  mode (genuinely similar, non-identical sequences winning on merit).
- **Reverted** to PR #6's cross-exclusion-only design for all four
  families; finding documented in each family's `families.yaml`
  description. The dedup fix itself is kept — it's a valid general
  safeguard for betL and any future `decoy_from_negatives` use, it just
  wasn't sufficient on its own here.
- No further action planned for proX/opuAC/opuBC/opuCC's remaining
  precision problem — it's now an open question without a known fix, not
  a pending decision.

### v7 full rebuild — verification result

Built and benchmarked all 43 families (`releases/v7`) against the last
clean `v6` build to verify the full sweep before committing.

- **Per family** (each weighted equally): DIAMOND 22 families improved
  (sum +0.817 F1) vs. 18 worsened (sum -0.413); HMM 25 improved
  (sum +0.706) vs. 17 worsened (sum -0.373). Net positive on both methods.
- **Read-volume-weighted** (pooling tp/fp/fn across all families, which
  better reflects real-world system performance): essentially **flat**.
  DIAMOND F1 0.823→0.821, HMM F1 0.684→0.684. DIAMOND caught 2,018 more
  true positives in absolute terms (real completeness gain, matching the
  sweep's actual goal) but picked up proportionally similar new false
  positives/negatives too.
- **Honest framing**: this is a completeness fix (captures real gene
  instances the panel was previously missing entirely), not a
  performance-optimization win — worth stating plainly rather than
  implying detection got "better" overall.
- **9 regressions above 0.03 F1** documented individually in
  `families.yaml` rather than chased: ectA (HMM), ktrA (DIAMOND), ktrD
  (HMM), mrpB (DIAMOND, -0.086, the largest), opuAB (HMM), opuCA (HMM),
  opuCB (both methods), trkA (DIAMOND, already known from Phase 2,
  confirmed at 0.920→0.887 in the full build vs. 0.920→0.879 in the
  earlier 6-family subset test).
- **Methodological finding**: a subset test isn't a fully reliable
  predictor of full-panel behavior — some families (ktrA, ktrD) showed
  different results in the full 43-family build than they did in the
  earlier 6-family subset test, because the full-panel build changes the
  combined DIAMOND/HMM database's composition beyond what a small subset
  can capture. Worth remembering for any future phase that relies on a
  subset test as a go/no-go signal.

### Phase 5 — trkH/ktrB/ktrD re-evaluation (PR #11)

Checked whether the numbered-paralog fix (Phase 2) or a decoy conversion
(Phase 4's approach) could address trkH/ktrB/ktrD's remaining negative-pool
contamination, per the plan's original Phase 5 framing.

- **Numbered-paralog fix didn't help**: purity-flag rates in v7 barely
  moved from v6 (trkH 27%→29%, ktrB 14.5%→15.5%, ktrD 2.4%→2.8% flagged)
  — confirms the same structural "unlabeled ortholog" pattern as Phase
  3's families, not a missing-symbol gap.
- **Decoy conversion ruled out**: trkH/ktrB/ktrD share PF02386 the same
  close-paralog way proX/opuAC/opuBC/opuCC shared PF04069 — exactly the
  clique shape Phase 4 showed catastrophically collapses DIAMOND recall.
- **No `pfam_model` fallback possible**: removing it from all three was
  the original fix for the byte-identical-HMM collision on this same
  trio (see PR #6's bug-fix history above).
- Both DIAMOND (0.265–0.293) and HMM (0.121–0.288) precision are low
  across all three, with no independent rescue via `profile_cascade.tsv`.
  ktrD's HMM is the worst in this cluster (0.121 precision, 0.904
  recall — calling almost everything positive), consistent with its
  already-noted small-n calibration risk (54 UniProt positives).
- Set `scope: annotate_only` on all three, same criteria as otsA/mrpC:
  removes the entire Trk/Ktr K+-uptake system from `osmotool profile`'s
  reported output (a bigger call than extending scope on 1-2 families,
  flagged and confirmed explicitly rather than applied by default) while
  keeping it built and searchable for `osmotool annotate` co-occurrence
  checks.
- Also backfilled a documentation gap: mrpC's `scope: annotate_only`
  change (from the prior commit) had been applied to `families.yaml` but
  never reflected in README's "Structural negative-pool contamination"
  section — fixed alongside trkH/ktrB/ktrD.

### Status at time of writing

- **Phase 1, 2, 3**: committed (`6d03487` on `fix-numbered-paralog-gaps`,
  on top of the original 6-family commit `69f80a4`), pushed, PR #9 and
  issue #8 updated with the final findings.
- **Phase 4**: disproven and reverted, documented in the same commit —
  proX/opuAC/opuBC/opuCC match PR #6's original cross-exclusion design.
- **Phase 5**: done — committed (`8b49d2d` on `phase5-trk-ktr-scope`),
  pushed, PR #11 open against `fix-numbered-paralog-gaps`.
- **Phase 6** (qc_scorecard-as-gate documentation for `pfam_model`): not
  yet started.

## This branch (`cluster-negatives-cdhit`) — negative-pool contamination & length-filter calibration

Branched off `main` after PR #13 merged. Two starting points: the branch's
namesake fix (CD-HIT wasn't deduplicating negatives, only positives, before
the train/test split) and a live question about whether the length-outlier
filter's median — computed from a *capped, sampled* negative fetch — was
itself part of why some families' negative pools perform poorly.

### CD-HIT clustering extended to negatives (`0038d82`)

`02_cluster_cdhit.sh` previously deduplicated only `*.positive.faa`.
Near-duplicate negatives splitting randomly across train/test lets a
`decoy_from_negatives` family's DIAMOND reference "see" a near-identical
sequence to one of its own benchmark negatives, inflating apparent decoy
performance — the same anti-leakage argument that already justified
clustering positives. Now clusters both labels for every family.

### Extra-positives merge step (`1e9c887`)

Added `01d_add_extra_positives.py`, wired in as build step 1d (between
length-outlier filtering and clustering): merges manually-curated hits from
a study's own assemblies (e.g. Bakta calls this repo's own detector missed
on mazG/mscS/trkA/trkH) into a family's positive set via
`extra_sequences/<family>.faa` (gitignored, local study data). Tested via
`families_bakta_test.yaml` against those four motivating families.

### Extra-positives merge + negatives-clustering: verified end-to-end (mazG/mscS/trkA/trkH)

Ran steps 1→1b→1c→1d→2 (fetch → purity filter → length filter → merge →
cluster) against `families_bakta_test.yaml`'s four motivating families in a
throwaway release dir, to confirm the merge step and the negatives-clustering
fix both behave as intended on real data before trusting them in a full
rebuild.

- Bakta-derived positives survive CD-HIT dedup at roughly a third to a half
  per family: mazG 162/455 (36%), mscS 111/254 (44%), trkA 279/821 (34%),
  trkH 119/253 (47%) — genuinely novel sequences at 90% identity, not just
  redundant restatements of what UniProt already had.
- Negative-set CD-HIT redundancy is much lower than positives' (1-9% vs.
  30-65%), consistent with a broad, taxonomically diverse hard-negative pool
  rather than one gene's own orthologs: mazG 140→138, mscS 615→590,
  trkA 551→499, trkH 694→667.

### cspA's 19% positive-length-outlier rate: real biology, not a fixable bug (diagnosis only)

Investigated why cspA's positive set (`pos_length=19%` in `qc_scorecard.tsv`)
flags far more than typical, since the panel-wide median-instability sweep
below doesn't cover the *positive* side or non-capped, `max_positive_override`
-sampled fetches like cspA's (5,000 of ~80,000 PF00313 members).

- 959/961 flagged sequences are too *long*, not fragments (flagged median
  154aa vs. kept median 68aa, up to 1,335aa). 86% (830/959) have no second
  Pfam domain annotated at all — not a fusion-partner case like mrpA/mrpB or
  otsA/otsB below — and the 129 that do have one show no common partner
  (DUF1294, ribosomal S30EA modulation, Excalibur Ca-binding, NYN, RNase,
  AhpC/TSA — a grab-bag). No dominant taxonomic skew either (Actinobacteria
  genera are ~20% of the no-second-domain group, spread thin otherwise).
- Conclusion: genuine, broadly-distributed long-form CSD-paralog biology
  that cspA's own `families.yaml` entry already says to capture ("expect
  multiple paralogs per genome... capture all matches") — the generic
  median-based length filter's default `max_ratio=1.5` is simply too tight
  for this family's real length heterogeneity, unlike the fusion
  contamination cases below, which are genuine data-quality problems.
- Calibrated what a fix would need: `max_ratio≈5.0` keeps 99.8% of cspA's
  positives (vs. 80.8% today), leaving only the 4 genuinely isolated extreme
  outliers (408-1,335aa) flagged. But `01c_check_length_outliers.py`'s
  `--max-ratio` is currently one global CLI flag applied to all 43 families,
  not a per-family override — acting on this would need a new
  `families.yaml` field (e.g. `length_max_ratio_override`), mirroring the
  `max_negative_override`/`max_positive_override` pattern already used
  elsewhere in this panel. **Not implemented** — diagnosis only.

### mrpA/mrpB: median instability disproven, real bug found instead (`5fd6e1f`)

Tested the original hypothesis — recalibrate a family's length-filter
median from a larger, uncapped fetch — on two documented bad performers.

- **mrpC**: negative-set median identical (102aa) whether fetched capped
  (n=1000) or fully uncapped (n=24,356). Not a length-filter problem — its
  low precision is the same "real orthologs UniProt never assigned a
  symbol to" structural contamination as Phase 3 documented for other
  families; ruled out as a length-filter fix target.
- **mrpB**: median *also* stable regardless of sample size (937–944aa) —
  but that stable number was itself wrong. Found ~68% of its negative pool
  was genuine mrpA+mrpB fused-ORF sequences under bare locus tags (same
  PF13244+PF20501+PF04039+PF00361+PF00662 architecture as mrpA's
  documented fusion cases, confirmed via direct UniProt lookup), invisible
  to the negative_query's gene-symbol exclusion and to 01c's
  positive-fetch-only fusion detection. Being the majority of the pool,
  they set the "normal" median, causing the filter to discard the true
  ~140aa hard negatives as outliers instead of the contamination.
- **Fix**: 01_fetch_refs.py now fetches Pfam domain evidence for the
  *negative* fetch too (previously positive-only) for `fusion_partner`
  families; 01c strips domain-confirmed fused-ORF sequences out of the
  negative pool *before* computing the median (length can't be used here
  — it's exactly what the contamination corrupts). 08d_build_fusion_refs.py
  picks these up as a bonus addition to the shared `mrpA_mrpB_fused`
  DIAMOND reference (they're genuine detection targets, not just
  decontaminated negatives). Confirmed: mrpB's negative median drops
  937aa→161aa; mrpA sees the same effect at smaller scale (~5%, since its
  negative anchor domain PF00361 is broader/more diluted). Regression
  -tested against all 84 other family/label combinations in the existing
  `v7` build — byte-identical.

### otsA/otsB: same bug found a second time (`fa6f229`)

Tested three more candidates (ktrA, otsB, opuCA) for median instability —
all three stable, ruling that out — but investigating *why* their medians
sat far from their own positive medians surfaced a second instance of the
mrpA/mrpB bug:

- **otsA/otsB**: ~41–59% of both families' negative pools are genuine
  bifunctional trehalose-phosphate synthase+phosphatase proteins
  (PF00982+PF02358 in one ORF, bare locus tags, some annotated
  "glucosylglycerol-phosphate synthase" — a closely related bifunctional
  enzyme with the same two-domain shape) — an undeclared analog of the
  mrpA/mrpB fusion, confirmed via direct UniProt lookup. Declared
  `fusion_partner`/`fusion_marker_pfam` for both, reusing the mrpA/mrpB
  mechanism with no code changes needed. Confirmed: otsA's negative median
  drops 724aa→487aa, otsB's drops 725aa→263aa (now matching its own
  positive median almost exactly).
- **ktrA**: dominant negative cluster is unrelated NhaP2/generic
  RCK-domain proteins sharing the broad PF02080 domain — genuine
  promiscuity, not a hidden fusion, no quick fix.
- **opuCA**: dominant cluster is `guaB`/IMP dehydrogenase, unrelated,
  sharing the broad CBS domain (PF00571) — same story, no quick fix.

### Pushing the fix upstream to query time (`4d8ca08`)

Negative queries for `fusion_partner` families now exclude the partner's
marker domain at fetch time (`NOT xref:pfam-<marker>`), not just post-hoc
— confirmed on mrpB: drops the raw population from 9,774 to 3,739, so the
same n=1000 fetch budget lands entirely on genuine hard negatives instead
of mostly-discarded contamination (previously only ~231/1000 survived).

Added a dedicated, uncapped combined-domain query per declared fusion pair
to comprehensively recover fusion candidates directly, replacing reliance
on incidentally spotting them. Getting the query right took two false
starts, both instructive:

- A naive "AND the two declared marker domains" breaks when a pair shares
  one marker value (mrpA/mrpB both declare PF13244) — substituting a
  family's own `negative_pfam` isn't safe either, since that domain can
  also be part of the *same* family's own standalone architecture
  (confirmed: `PF00361 AND PF13244` just re-matched ~8k ordinary mrpA
  orthologs, not fusions).
- Fixed by verifying each candidate domain empirically against each
  family's own already-fetched `positive.domains.tsv` (self-intrinsic to
  one side, foreign to the other) before combining them. Even that can't
  resolve every pair: mrpA is itself a large multi-domain protein (own
  median ~800aa) whose size range overlaps the fused-with-mrpB range
  (~940aa), so no domain pair can be verified safe for it — the dedicated
  fetch correctly skips that pair in favor of the already-validated
  incidental-discovery + post-hoc pre-filter. otsA/otsB's dedicated fetch
  works cleanly (1,390 candidates via `PF00982 AND PF02358`).
- Also fixed two stale-file bugs found while testing this (in both
  `01_fetch_refs.py` and `08d_build_fusion_refs.py`): a re-run that
  decides to skip a pair now deletes any dedicated-fetch/merged-fusion
  output left over from a previous (or buggy) run, instead of silently
  leaving it for the next step to trust.

### Full-panel sweep: mazG is the only genuine median-instability case (`c752e7e`)

Systematically checked every family whose negative fetch hits the default
1000-sequence cap (38 of 43) for capped-vs-larger-fetch median instability
— the original hypothesis this branch started from.

- **Result: mazG is the only case in the entire panel.** Negative-set
  median is 215aa at the default n=1000 cap vs. 116aa fetched at n=3000 —
  a real, sample-size-driven effect. Every other family tested moved by
  ~2% or less between those two fetch sizes (largest: cspA at 2.5%).
- Two families hit a transient UniProt read-timeout during the sweep
  (opuAB, opuBB) — re-fetched individually, both confirmed stable (~1%
  change) once recovered.
- **Fix, two parts** (fixing the median alone wasn't enough — see below):
  new `max_negative_override` field (mirrors the existing
  `max_positive_override` pattern) raises mazG's own fetch cap to 3000, so
  a normal build actually pulls enough sequences near the true center;
  confirmed that applying the correct 116aa median to the *old* capped
  1000-sequence fetch left only 5 survivable negatives, since that small
  sample barely contained anything near the true center in the first
  place. New `negative_median_override` field pins the value itself (01c
  reads it instead of recomputing), as protection against the value
  drifting if UniProt's underlying data changes over time. Together: 689
  well-calibrated negatives survive now, vs. 140 (wrong median) or 5
  (right median, too-small sample) before.
- Given only one family out of the whole panel needed this, it's
  implemented as a documented one-off (two families.yaml fields on mazG's
  entry), not a general mechanism.

### Screened mrpF/mrpG/gshA/gshB/trkH/ktrB/ktrD for the same fusion bug — none found

Checked whether any of the panel's other documented "structural
negative-pool contamination" families were actually an undeclared fusion
pair like otsA/otsB, using data already on hand (no new fetches needed):
compared each family's own negative-set median against its positive-set
median. The confirmed fusion cases showed the negative median dramatically
larger than the positive one (mrpB ~6.6x, otsB ~2.7x); all seven candidates
here sit within 0.92x–1.06x — mrpF 93/93aa, mrpG 119/121aa, gshA 481/521aa,
gshB 334/317aa, trkH 462/483aa, ktrB 482/454aa, ktrD 482/453aa. None show
the signature. Their documented contamination is most likely the plain
"real orthologs UniProt never assigned a symbol to" pattern from Phase 3,
with no quick fix available — consistent with mrpC's finding above.

### Environment note

Running a real build+benchmark (rather than just the fetch/filter steps)
needed packages missing from the `osmotool` conda environment relative to
`environment.yml`: `trimal`, `wgsim`, `biopython`, `scikit-learn`,
`insilicoseq`. Installed via conda-forge/bioconda into that environment
(all available there; no `pip`/PyPI access was needed). Separately, `mafft`
needs real system temp-directory access that this session's sandboxed
shell blocks by default — commands invoking it need the sandbox disabled
for that step specifically.

### Full v8 rebuild — benchmark validation, a real bug, and a partial revert (`16bf1d5`, `0151ed6`)

Built and benchmarked all 43 families (`releases/v8`) to confirm the
mrpA/mrpB/otsA/otsB/mazG fixes actually move DIAMOND F1/precision, not
just the calibration median as confirmed up to this point.

- **Found and fixed a real, pre-existing bug along the way, unrelated to
  this branch's own work**: the build itself completed cleanly through
  all 43 families, but the benchmark step failed — `osmotool profile`'s
  `DATABASE` argument is an unpacked release *directory* (it finds
  `osmo_refdb.dmnd`, `hmms/osmo_refdb.hmm`, etc. inside it by fixed
  names), while `run_pipeline.sh`/`10_run_benchmark.sh` were passing the
  `.dmnd` *file* path directly — an interface drift between this repo and
  `osmotool` that would have broken every future benchmark run. Fixed by
  passing the release directory itself; confirmed with a single-sample
  test before resuming the full benchmark.
- **Panel-wide (read-volume-weighted)**: essentially flat, as expected
  since only 5 of 43 families changed — DIAMOND F1 0.821→0.820, HMM
  0.684→0.690.
- **Per family**, the picture was mixed, not a clean win:
  - **mazG**: clear win, exactly as targeted — DIAMOND F1 0.855→0.883
    (precision 0.873→0.992, FP 1019→56), HMM F1 0.822→0.854.
  - **otsA**: real improvement — DIAMOND F1 0.655→0.700 (precision
    0.518→0.576), HMM F1 0.579→0.612.
  - **mrpA**: essentially flat — DIAMOND F1 0.714→0.701.
  - **otsB**: regressed — DIAMOND F1 0.829→0.773 (precision
    0.902→0.761), HMM precision 0.839→0.711.
  - **mrpB**: collapsed — DIAMOND F1 0.587→0.272 (precision
    0.888→**0.187**), recall actually rose slightly (0.438→0.498).
- **Root-caused the two regressions, not just reverted them**:
  `08b_calibrate_diamond_cutoffs.py` sets its cutoff from the negative
  set's own score distribution. The fusion-contamination sequences
  removed by this branch's fix, despite being conceptually wrong as
  ground truth, were genuinely domain-similar to these families'
  positives (real shared domain content, just fused to the partner's),
  so they forced a strict calibration threshold. Stripping them out left
  only very dissimilar "easy" negatives that don't constrain the cutoff
  enough — it drifted permissive and let through more false positives
  from elsewhere in the 43-family DIAMOND db (recall rising alongside a
  precision collapse is the signature of a loosened threshold, not a
  worse reference set). Same mechanism confirmed independently on both
  otsB and mrpB, at proportional severity to each one's contamination
  fraction (otsB ~59%, milder regression; mrpB ~68%, the collapse).
- **Partial revert**: removed `fusion_partner`/`fusion_marker_pfam` from
  mrpB and otsB specifically — back to their original `negative_query`,
  no negative-side pre-filter or query-time exclusion. mrpA/otsA's fixes
  are kept (confirmed via live fetch that each family's own
  decontamination runs independently of its partner's declaration, so
  reverting mrpB/otsB doesn't touch mrpA/otsA and vice versa). Accepted
  side effect on mrpB: this also disables its own positive-side
  oversized-candidate check (the code ties both to the same
  declaration) — mrpA's own positive-side check and taxonomic evidence
  are unaffected and still catch most such candidates.

### Revert confirmed on a 4-family subset rebuild; a real mrpA puzzle explained (`a6da7c6`)

Re-ran build+benchmark for just mrpA/mrpB/otsA/otsB (`releases/revert-confirm`)
to confirm the partial revert actually recovers mrpB/otsB rather than
trusting the code-level revert alone.

- **mrpB and otsB fully recovered, both edging past their original v7
  baseline**: mrpB DIAMOND precision 0.187→0.857 (F1 0.272→0.618, vs. v7's
  0.587), otsB precision 0.761→0.936 (F1 0.773→0.843, vs. v7's 0.829).
- **mrpA looked worse in this subset (DIAMOND F1 0.701→0.639) despite its
  own fix being unchanged — root-caused, not left as a mystery**: 3,330 of
  mrpB's 8,600 reverted-negative-set reads (38.7%) score as genuine `mrpA`
  hits, confirmed directly from `mrpB.negative.gene_counts.tsv`. This is
  mechanically correct, not a DIAMOND error: reverting mrpB restored the
  same mrpA+mrpB fusion contamination found earlier in this branch, and
  those sequences really do contain mrpA's own domain content. The
  benchmark's ground truth labels them "mrpB negative" (so a `mrpA` call
  counts as a false positive there), even though the call is sequence-
  content-correct. A side effect of reverting mrpB, not a flaw in mrpA's
  own kept fix.
- Otherwise a clean, apples-to-apples confirmation that each family's own
  decontamination runs independently of its partner's declaration, exactly
  as designed.

### Read-simulation reproducibility fix (`a6da7c6`)

While comparing runs, noticed the same family's total simulated-read count
wasn't even identical between two supposedly-comparable builds. Root
cause: neither `wgsim` (`-S`) nor InSilicoSeq (`--seed`) were ever passed
an explicit seed, so — unlike every other random draw in this pipeline —
the actual simulated reads (counts, sequencing errors, sampled positions)
differed between runs even against byte-identical held-out input
sequences. Fixed: each family/label's simulator call now gets a
deterministic seed drawn from `09_simulate_reads.py`'s own seeded rng.
Confirmed via a clean two-run comparison: byte-identical FASTQ output.
Live UniProt fetches remain a separate, unaddressed source of run-to-run
variation (accession sampling is seeded, but the underlying database can
change between runs done at different times).

### Status at time of writing

### Full `v8` rebuild redone with the partial revert + seed fix in place — final confirmation

The first `v8` rebuild (above) predated the mrpB/otsB partial revert and
the read-simulation seeding fix. Re-ran the complete 43-family build+
benchmark from scratch with both in place, to get one clean, final,
reproducible number for the whole branch's work rather than reasoning
from a subset test and a stale full-panel run.

- **Panel-wide (read-volume-weighted)**: genuine net improvement now,
  not the earlier flat/mixed result — DIAMOND F1 0.821→0.825, HMM
  0.684→0.693.
- **mazG, otsA**: unchanged clean wins (DIAMOND F1 0.855→0.879,
  0.655→0.699).
- **otsB**: fully recovered *and* edges past its v7 baseline in the full
  panel too (DIAMOND F1 0.829→**0.840**) — confirms the revert isn't a
  subset-only effect.
- **mrpA**: actually improved here (DIAMOND F1 0.714→**0.748**, better
  than even the original buggy-fix `v8`'s 0.701) — confirms the
  subset-test regression really was the subset-composition artifact
  diagnosed above (mrpB's contamination dominating a tiny 4-family read
  pool), not a real problem with mrpA's own kept fix at full-panel scale.
- **mrpB**: essentially back to its v7 baseline (DIAMOND F1 0.587→0.569,
  a 0.018 difference) — within normal fetch-to-fetch noise from live
  UniProt data drifting slightly between runs, not a residual issue. The
  revert did what it was supposed to.

This is the number to cite for this branch's net effect: five families
touched, three genuine wins (mazG, otsA, otsB), one improvement that only
showed up at full-panel scale (mrpA), one full recovery to baseline
(mrpB), and a real panel-wide improvement rather than the flat trade-off
the first attempt suggested.

### Status at time of writing

- CD-HIT negatives-clustering fix, extra-positives merge step, the
  mazG fix, the benchmark-path bug fix, the mrpA/otsA fixes (kept), the
  mrpB/otsB partial revert, and the read-simulation seeding fix are all
  **committed** on `cluster-negatives-cdhit`, not yet pushed or opened as
  a PR.
- **Full-panel benchmark validation is done** — the final `v8` rebuild
  (immediately above) confirms every fix and the revert on real 43-family
  numbers, not just the subset test or the pre-revert build. Nothing
  outstanding on the negative-pool contamination work itself.
- Extra-positives merge + negatives-clustering verified end-to-end against
  the four motivating families in a throwaway subset build (see above) —
  confirms the mechanism works on real data, not yet folded into a full
  43-family rebuild (the final `v8` rebuild used `families.yaml` directly,
  without an `extra_sequences/` directory present).
- cspA's positive-length-outlier mismatch is diagnosed and a fix calibrated
  (`max_ratio≈5.0`), but **not implemented** — needs a new per-family
  override field in `01c_check_length_outliers.py`/`families.yaml` that
  doesn't exist yet.
- `docs/FAMILIES_SUMMARY.md` is untracked in the working tree, predates
  this session's work, no action taken.
- Not yet done: pushing `cluster-negatives-cdhit` and opening a PR.

## This branch (`add-refseq-positives`) — RefSeq as a second positive source

Follow-up to the extra-positives merge above. Benchmarking the merged `v8`
build showed the Bakta-study merge alone doesn't close the recall gap it
targeted: reads simulated from the held-out (test-split) study-derived
sequences recall far worse than reads from ordinary UniProt sequences (e.g.
mazG DIAMOND 23.2% vs. 80.7%, trkH DIAMOND 6.5% vs. 66.2%) -- the merged
sequences are genuinely divergent from what the reference already had (the
same reason they survived CD-HIT clustering as non-redundant in the first
place), and divergence is exactly what a single calibrated score threshold
struggles to generalize to.

Researched how Bakta itself annotates these families successfully, for
comparison: its DIAMOND search runs at a much more permissive 80%
coverage/50% identity bar (vs. one calibrated per-family cutoff here), and
its reference corpus (UniRef90/UniRef50/IPS, built from clustering all of
UniProt+RefSeq) is far broader than any single gene-symbol-anchored UniProt
query. That pointed at RefSeq as a way to broaden the *positive* corpus --
a different lever than the already-disproven "fetch more negatives"
(`max_negative_override` on mrpF/mrpG, Phase 1) -- worth testing directly
rather than assuming.

### Feasibility check across 15 families

Confirmed mazG/mscS/trkA/trkH's UniProt positive populations are already
fully exhausted (fetched count == UniProt's live total for all four --
no fetch cap left to raise), so growing them further needs a source outside
UniProt's own gene-symbol curation. Ran a CD-HIT-2D novelty check (RefSeq
sample vs. the existing UniProt+study-merged pool, 90% identity) plus a
manual product-description spot-check, across those four and 11 more
families flagged earlier in this branch's negative-survival ranking:

- **9 families confirmed clean and worth adding**: mazG (41.0% novel, RefSeq
  population 45,441 vs. UniProt's 5,340), mscS (26.4%), trkA (25.2%), murB
  (32.8%, 74,987 vs. 12,332), mscL (42.8%), otsB (65.6%, largest novelty
  rate in the panel), otsA (29.8%), gshB (18.4%), mrpG (76.6%). All showed
  ~94-100% correct product-description annotation on the novel subset --
  no cross-contamination.
- **Skipped for low marginal value (clean but not worth it)**: trkH (4.4%
  novel), ktrB (12.8%, tiny 415-member RefSeq population to begin with),
  ktrD (0% novel -- RefSeq offers nothing this pool doesn't already have).
- **Skipped for gene-symbol collision (real contamination found, not
  hypothetical)**: mrpC -- only 32% of its "novel" RefSeq hits are
  genuinely mrpC (Na+/H+ antiporter subunit C); 45% are "MR/P fimbria
  usher protein MrpC" (a *Proteus* fimbrial gene) and 23% are "Crp/Fnr
  family transcriptional regulator MrpC" (a *Myxococcus* transcription
  factor) -- the bare symbol "MrpC" is shared by three unrelated gene
  families. mrpF -- a smaller-scale version via its own documented PhaF
  alias (PR #7), which also names an unrelated polyhydroxyalkanoate-granule
  protein in RefSeq. mrpB -- "DUF1883 domain-containing protein MrpB"
  confirmed (direct lookup) to be an unrelated small protein (DUF1883 +
  PPC + bacterial SH3-like domains, found in e.g. *C. difficile*), not a
  real mrpB variant. None of the three siblings extended past mrpG.

### Implementation (`pipeline/01e_add_refseq_positives.py`)

New opt-in `families.yaml` field, `refseq_gene_symbols`, set only on the 9
confirmed-clean families -- deliberately not a general mechanism, given the
collision risk found above. New build step 1e (between the extra-positives
merge and CD-HIT clustering) fetches NCBI RefSeq bacterial proteins tagged
with those gene symbols and merges them into the fetched positive set:

- Fetches the full matching UID list first (lightweight -- UIDs only) and
  draws a reproducible random sample from it, rather than trusting
  esearch's own non-randomized default order -- the same class of bias
  01_fetch_refs.py already had to fix for UniProt (issue #4).
- Skips RefSeq entries explicitly marked ", partial" at fetch time (a free
  fragment signal); everything else still goes through the normal
  length-outlier filter downstream.
- Retry-with-backoff on NCBI's eutils calls -- confirmed necessary, not
  speculative: hit a real transient 502 and two "response ended
  prematurely" errors during testing, all recovered cleanly on retry.
- Idempotent re-run (dedup by RefSeq accession), matching 01d's pattern.
- Verified end-to-end in a throwaway subset (`releases/refseq-test`,
  removed after): 7,350 RefSeq sequences merged across all 9 families
  (1,435-1,492 sampled per family from each family's ~1,500-target cap),
  clean sequence content (standard 20-aa alphabet, sane length ranges), no
  errors after the retry fix.

### Full build+benchmark validation, one family at a time (mrpG, otsA, mscS, otsB)

Ran single-family build+benchmark subsets (`make_family_subset.py`) with the
RefSeq merge active, each compared directly against its own `v8` baseline
(same code otherwise, since `v8` predates `01e` entirely) -- the actual test
the sequence-merge-level validation above couldn't answer on its own.

- **mrpG** (chosen as the worst v8 performer, F1 0.541/0.442, and the
  highest RefSeq novelty rate of the 9 at 76.6%): DIAMOND F1 0.541->0.465
  (worse), HMM F1 0.442->0.472 (barely better). Splitting recall by origin
  showed why this one didn't move much either way: RefSeq-origin recall
  wasn't dramatically worse than the rest for either method (DIAMOND 32.7%
  vs 39.3%, HMM 38.2% vs 42.0%) -- mrpG's badly-performing status turned out
  to be its already-documented `pfam_ga_review_needed` negative-pool
  calibration problem (Phase 3), not a positive-pool-narrowness problem, so
  growing the positive pool was the wrong lever for this specific family.
- **otsA** (next-worst cap-limited family without mrpG's calibration flag):
  DIAMOND F1 0.699->0.631 (worse), HMM F1 0.607->0.657 (genuinely better --
  both precision and recall improved). RefSeq-origin recall: DIAMOND 54.2%
  vs 77.8% (real gap), HMM 87.3% vs 86.5% (no gap, explains the clean HMM
  win).
- **mscS**: the one exception -- DIAMOND F1 0.495->0.512 (slightly better),
  HMM F1 0.455->0.505 (clearly better). Even here, HMM's gain came from a
  real recall edge on RefSeq-origin reads (77.8% vs 64.1%) that DIAMOND
  didn't share (57.9% vs 65.9%, still a deficit) -- DIAMOND's improvement
  came despite the dilution, not because RefSeq-origin sequences were
  suddenly easy for it too.
- **otsB** (highest RefSeq novelty rate, 65.6%): DIAMOND F1 0.840->0.717
  (worst regression yet, -0.123), HMM F1 0.662->0.668 (flat/noise).

**Consistent pattern across all 4**: HMM's F1 improved in every case (3
meaningfully, 1 negligibly); DIAMOND's F1 got worse in 3 of 4, including the
largest single regression seen anywhere in this branch's testing (otsB,
-0.123). The mechanism matches the Bakta-merge finding and the earlier
Bakta-methodology research: profile-based matching (HMM) tolerates the
added divergent sequences; identity-based best-hit search (DIAMOND) mostly
doesn't.

### Per-method reference split: full merge for HMM, UniProt-only for DIAMOND (`08_build_diamond_db.sh`)

Direct response to the pattern above, proposed and confirmed rather than
assumed: since HMM and DIAMOND consistently "need different things," build
each from a different input instead of forcing one merged pool to serve
both.

- `08_build_diamond_db.sh` now excludes `_study`/`_refseq`-tagged sequences
  from the DIAMOND reference specifically -- `positive.train.faa` itself is
  untouched (04_align_trim.sh already consumed the full file for HMM's
  alignment before this step runs), so nothing upstream needs to change.
  Extended to also exclude `01d`'s `_study` tag, not just `01e`'s `_refseq`
  one -- the same mechanism applies there (confirmed earlier: Bakta-origin
  recall for DIAMOND was even worse than RefSeq-origin's, e.g. mazG 23.2%
  vs 80.7%), even though the request that prompted this was RefSeq-specific.
- Negatives are untouched: there's no RefSeq/study-equivalent negative
  source (`01d`/`01e` only ever touch positives), so both methods still
  share the same UniProt-only negative pool -- nothing to split there yet.
- Implemented as a small `awk` filter (drop whole records by tag, keep
  everything else) rather than a new pipeline step, since `03`/`04` already
  read the unfiltered file correctly and only DIAMOND's own db-building
  step needed to diverge. Verified directly against real `v8` data
  (mazG's 2,629-sequence train set, 134 `_study`-tagged): filtered output
  is exactly the 2,495 non-`_study` records, byte-identical to the
  originals, nothing else disturbed.

### Full 43-family `v9` rebuild -- panel-wide validation

Built and benchmarked all 43 families (`releases/v9`) against `v8`, the same
"validate the whole branch before merging" discipline this project used for
its own `v7`/`v8` rebuilds -- the single-family tests above show the
mechanism works, but only a full-panel build confirms it at the scale this
will actually ship at.

- **Panel-wide (read-volume-weighted)**: essentially flat, as expected --
  DIAMOND F1 0.825->0.827, HMM F1 0.693->0.693. This branch's changes only
  touch 9 of 43 families, so a flat panel-wide number isn't a null result,
  it's the correct outcome when most of the panel is untouched.
- **The 9 RefSeq-enabled families, isolated** (the actual test): DIAMOND F1
  0.854->0.858 (+0.003), HMM F1 0.682->0.683 (+0.001), both combined
  read-volume-weighted across just these 9. Per-family, DIAMOND improved or
  stayed flat in all 9 (mazG +0.014, mscS +0.035, trkA +0.050, murB -0.001,
  mscL +0.009, otsB -0.003, otsA +0.003, gshB -0.001, mrpG +0.030) -- no
  regressions at all, a real change from the *unsplit* single-family tests
  above where otsA/otsB/mrpG all regressed (up to -0.123 for otsB). That
  regression recovery is the actual confirmation the DIAMOND/HMM split
  works, not just a plausible mechanism. HMM stayed mostly positive too (6
  of 9 improved, murB -0.019 the largest dip, otherwise noise-level).
- **9 other (non-RefSeq-enabled) families moved by more than 0.03 F1** in
  one direction or the other on one or both methods (ktrD, opuCB, opuBA,
  mrpC, opuCC, mrpB, ktrA, opuAA, opuCA, opuBB) -- none of these declare
  `refseq_gene_symbols` or have an `extra_sequences/` file, so nothing in
  this branch touches them directly. Consistent with this project's own
  documented finding (Phase 2's numbered-paralog sweep, `v7` rebuild
  notes): a full-panel build shifts the combined DIAMOND/HMM database's
  composition for every family, not just the ones a change directly
  targets, plus ordinary live-UniProt-data drift between fetches done on
  different days. Not chased further -- same "document, don't chase"
  discipline used throughout this project's history for unexplained
  full-panel drift below the scale of a real regression.

### Reproducibility audit: broken provenance, environment drift, and a real murB fix found along the way

Prompted by "this feels a bit all over the place" -- a fair read of a
branch that had accumulated a live-UniProt fetch, a live-NCBI fetch, an
ad-hoc hand-patched conda env, and several single-family throwaway builds.
Checked what was actually broken rather than just reorganizing.

- **UniProt release tracking was silently broken.** `get_uniprot_release()`
  hit `rest.uniprot.org/utils/release`, which 404s (confirmed directly) --
  every `manifest.tsv` row across `v7`/`v8`/`v9` recorded `unknown` instead
  of the real release, defeating the field's whole purpose. Fixed to read
  the `X-UniProt-Release` response header off an ordinary search request
  instead (present even on the cheapest possible query, `size=0`) --
  verified live: returns `2026_02`.
- **RefSeq fetches had no provenance at all.** `01e`'s manifest didn't even
  have a fetch-date column. Added `date_fetched`, the practical equivalent
  of UniProt's release field given NCBI has no discrete RefSeq release
  number.
- **The environment used for every benchmark this branch ran was not the
  one this repo documents.** Everything (mrpG/otsA/mscS/otsB, the full
  `v9` rebuild, the murB test below) ran through a pre-existing, hand-patched
  `osmotool` conda env, not a clean build from `environment.yml`. Built and
  verified the actual Docker image for the first time this session:
  `mafft` 7.526, `trimAl` 1.5.rev1, HMMER 3.4, DIAMOND 2.2.4, CD-HIT 4.8.1,
  `wgsim`, full Python stack, and `osmotool` all present and working --
  `pip show osmotool` confirms the exact pinned commit (`0.1.dev3+g4cbc3ba8d`
  matches the Dockerfile's `@4cbc3ba`). **Found a real discrepancy in the
  process**: the ad-hoc env's `osmotool` was a *different* commit
  (`0.3.1.dev2+gf1465475d.d20260727`, not `4cbc3ba`) -- every benchmark
  number produced on this branch so far used an unpinned, drifted build,
  not the one the project actually documents and pins for reproducibility.
  Numbers are still internally consistent (compared against each other and
  against `v8`, built the same way), but a future full rebuild through
  Docker is the one to treat as authoritative, not `v9` as currently built.
- **Added `.dockerignore`** (didn't exist before): the Docker build context
  was transferring all of `releases/` -- 11GB down to ~16MB, unrelated to
  correctness but a real drag on ever actually using Docker day-to-day.

### murB: fixed for real, found while investigating the reproducibility audit above

Re-examining murB's `insufficient_negative_data` status (the worst in the
43-family panel, see the survival-rate ranking earlier in this branch)
surfaced that its negative anchor, PF02873 (MurB_C, the catalytic
C-terminal domain), is *too* MurB-specific to serve as a negative pool at
all -- 41 domain architectures, essentially gene-specific -- so
`xref:pfam-PF02873 NOT gene:murB` wasn't selecting confusable paralogs, it
was selecting real MurB orthologs UniProt never got around to curating a
gene symbol for (confirmed: 53 total UniProt candidates, 75% flagged by
01b as likely mislabeled true positives, 13 left standing, hence the
`insufficient_negative_data` flag). The same narrow-domain pattern used
elsewhere in this file (proP, gshB) doesn't fail this way because those
genes have real distinct paralogs sharing their narrow domain; murB, a
near-universal single-function housekeeping enzyme, doesn't.

- **Tested rather than assumed**: switched `negative_query`/`negative_pfam`
  to PF01565 (FAD_binding_4, the N-terminal domain MurB shares with a
  broad class of unrelated FAD-linked oxidoreductases) and ran a full
  single-family build+benchmark against the `v9` baseline.
- **Calibration**: `insufficient_negative_data` (13 negatives, no workable
  cutoff) -> `clean_separation`, F1=1.000, cutoff=34.70, 61 usable
  negatives.
- **Real detection improvement, not just a calibration-status fix**:
  DIAMOND F1 0.936->0.941, HMM F1 0.705->**0.736** (recall 0.545->0.582,
  precision unchanged at ~1.0 for both methods -- a genuine gain on both
  axes, not a tradeoff).
- Kept anyway despite PF01565 being the "broad, promiscuous fold" pattern
  this file usually avoids (ectA/ectB/galE): every real MurB ortholog still
  carries PF02873 too, so 01b's purity filter still catches true positives
  pulled in under the broader anchor -- the promiscuity that would be a
  contamination risk elsewhere is exactly what supplies murB with negatives
  it never had access to before.

### osmotool Dockerfile pin was the stale side, not this branch's ad-hoc env

Follow-up on the "open item" above: checked the actual osmotool git history
(local checkout one directory up from this repo) instead of just flagging
the discrepancy and moving on.

- The Dockerfile's pin (`4cbc3ba`) was ~25 commits behind osmotool's
  current `main`. The ad-hoc env used for every benchmark on this branch
  (`f1465475d`) was only 2 commits behind current HEAD -- i.e. it was the
  Dockerfile pin that had drifted stale, not this branch's testing
  environment.
- One of those 25 commits (`c097a43`) fixed `quantifier.KNOWN_FAMILIES`
  from a hardcoded ~15-family list to the full v7 43-family panel --
  **mrpG, gshB, and trkA (this branch's own RefSeq/split test subjects)
  weren't in the old pin's known-family list at all.** A Docker build from
  the pre-bump Dockerfile would not have handled several of the families
  this branch actually validated.
- Bumped the Dockerfile's pin to `2826ebc` (current osmotool `main`,
  confirmed matching `origin/main`) and rebuilt. Verified: `osmotool
  --version` now correctly reports `0.4.0` (previously fell back to "see
  pyproject.toml" -- also fixed upstream in this commit range, not
  something broken in this repo), and the CLI's own `--help` no longer
  hardcodes the stale 15-family description. Rebuild was fast (~12s vs.
  the initial ~500s) confirming `.dockerignore` works as intended on a
  cache-hit build.

Net correction to the earlier framing: this branch's benchmark numbers
(mrpG/otsA/mscS/otsB, the full `v9` rebuild, murB's PF01565 test) were
already run against an environment close to current osmotool, not a
drifted/wrong one -- the Dockerfile just hadn't caught up yet. No need to
re-run any of this branch's validation; the Docker image is now the one
that matches what was actually tested, not the other way around.

### mrpE: catastrophic F1 root-caused as a self-defeating benchmark, not a detection problem

Investigated as the worst-performing family still in default `profile` scope
(`v9` `results/metrics/summary.tsv`: DIAMOND F1 0.262, HMM F1 0.196,
both driven by huge false-positive counts despite `qc_scorecard.tsv` showing
`hmm_status=pfam_ga_review_needed`, the same "structural, deferred to
Pfam's GA cutoff, considered fine" status already validated for mrpB/gshB).
That framing turned out to only justify the calibration-status flag, not
the real detection performance.

- **Root cause, confirmed live against UniProt, not assumed**: PF01899
  ("MNHE" in Pfam's own naming) has no Pfam clan or InterPro overlap
  (`set_info`/`overlaps_with` both null via the InterPro API) -- unlike
  murB's PF02873, there's no broader shared fold to re-anchor negatives on.
  A 500-entry live sample of the old negative pool (PF01899 minus
  gene:mrpE*, 9,839 total members vs. only 569 with a curated mrpE/mnhE
  symbol) showed ~99% are genuine MrpE/MnhE orthologs under a bare locus
  tag ("Na+/H+ antiporter subunit E", "Cation:proton antiporter", etc.) --
  the domain is essentially diagnostic for the gene, so the standard
  "same domain, no gene symbol" negative-pool recipe used everywhere else
  in this panel selects mislabeled true positives here, not hard
  negatives.
- **Confirmed via isolated single-family build+benchmark**
  (`make_family_subset.py`, same code both sides): swapping to a random,
  domain-unrestricted bacterial background (`taxonomy_id:2 NOT
  xref:pfam-PF01899`) dropped false positives to **zero** for both DIAMOND
  and HMM (precision 1.0/1.0) -- mrpE was never hard to detect; the old
  benchmark was measuring the negative pool's own label noise, not the
  model.
- **Recall gap addressed too**: removing `pfam_model: PF01899` (custom-
  built HMM, locally F1-calibrated instead of trusting Pfam's raw
  full-length-tuned GA cutoff) recovered real recall on top of the
  precision fix -- mrpE is a short protein (~150-250aa), so short
  simulated reads only cover a fragment of the domain and often legitimately
  scored below Pfam's 27.0-bit GA line. Local calibration (157.2 bits on
  whole-protein scoring) recovers more of them.
- **Combined result**: DIAMOND F1 0.239->0.75, HMM F1 0.19->0.704 (both
  isolated-build numbers, directly comparable pre/post since only
  `negative_query`/`pfam_model` changed). Adopted in `families.yaml`; see
  mrpE's entry for the full writeup, including the one known gap (PhaE,
  the sole real close paralog, is itself Pfam-integrated into PF01899 so
  it's excluded from the new negative pool same as everything else under
  the domain -- untested against that one specific confusable case).
- **Same mechanism confirmed present in mrpB and mrpF** via live UniProt
  sampling of their own negative pools (PF04039, PF04066): ~9,000-9,800
  candidates each, virtually all genuine MrpB/MrpF orthologs under a bare
  locus tag ("Monovalent cation/H+ antiporter subunit A/B/F", etc.), no
  real diversity. mrpG wasn't separately sampled but shares the same
  `pfam_model` + `pfam_ga_review_needed` shape. gshB, by contrast, looks
  genuinely different on inspection: its negative pool is small (246
  candidates, matching `qc_scorecard.tsv`) but shows real minority
  diversity (ATP-grasp domain proteins, a RimK-like ribosomal-protein-
  modification enzyme) alongside the majority mislabeled-gshB
  contamination -- not the same "domain is essentially the gene"
  situation, not redesigned. See the mrpB/mrpF/mrpG section below for how
  those three were actually fixed and validated.

### mrpB, mrpF, mrpG: same redesign as mrpE, validated as a group instead of in isolation

mrpE's isolated single-family test methodology (above) deliberately
controls for the "cross-family competition changes the numbers" confound
found earlier on ectB -- but it has its own blind spot, and mrpB's own
history already proved it's not hypothetical: mrpB's `families.yaml` entry
documents an earlier fix (removing fusion-contaminated sequences from its
negative pool, a much smaller decontamination than the redesign here) that
was tested via a full `v8` panel rebuild and caused DIAMOND precision to
collapse (0.888->0.187) -- `08b_calibrate_diamond_cutoffs.py` sets its
cutoff from the negative set's own score distribution, so "easier"
negatives calibrated a more permissive threshold that then admitted false
positives from *other* families in the full panel. An isolated
single-family test cannot detect that mechanism at all, since there's no
other family present for a false positive to come from -- which meant
mrpE's own fix (already adopted above) was validated against the wrong
kind of test for this specific risk.

- **Test design**: a 7-family subset (`make_family_subset.py --only
  mrpA,mrpB,mrpC,mrpD,mrpE,mrpF,mrpG`), built twice -- once with
  everything as currently committed (mrpE already fixed, mrpB/mrpF/mrpG
  still original), once with mrpB/mrpF/mrpG *also* redesigned the same
  way as mrpE (random background negative pool, `pfam_model` removed).
  mrpA/mrpC/mrpD kept unchanged in both builds specifically to provide
  real, biologically-related (same Mrp/Mnh complex) competing sequences
  in the DIAMOND db -- a sharper, cheaper probe for the specific
  calibration-drift risk than either total isolation or a full 90-minute
  43-family rebuild.
- **Result: the calibration-drift risk did not recur.** mrpB/mrpF/mrpG
  all went to exactly 0 false positives on both methods (mrpB DIAMOND
  19->0, HMM 718->0; mrpF DIAMOND 704->0, HMM 819->0; mrpG DIAMOND
  765->0, HMM 745->0), with large genuine F1 gains across all three
  (mrpB 0.649->0.675 DIAMOND / 0.248->0.692 HMM; mrpF 0.370->0.718 /
  0.273->0.743; mrpG 0.611->0.712 / 0.496->0.614). The unchanged sibling
  families showed no false-positive inflation: mrpA's numbers were
  identical between the two runs, mrpC/mrpD moved by under 0.04 F1 --
  matching this project's own documented "ordinary full-panel
  database-composition noise" (Phase 2 numbered-paralog sweep, `v7`
  notes), not the >0.4-F1 collapse pattern the earlier fusion-fix
  regression produced. Plausible reason the risk didn't recur here:
  mrpA/C/D are different Pfam domains entirely (PF00361, PF00420 vs.
  PF04039/PF04066/PF03334) with no real sequence homology to
  mrpB/F/G for a loosened cutoff to accidentally admit, unlike the
  fusion case, where the removed content was genuine shared PF04039
  sequence that had been anchoring the old threshold.
- **mrpE re-validated as a side effect**: this same 7-family test also
  exercises mrpE (previously isolation-only), with real mrpA/mrpC/mrpD
  competition present. Held up cleanly (F1 0.708->0.716 DIAMOND,
  0.69->0.727 HMM across the two 7-family runs) -- the isolation-only
  validation gap flagged above is now closed.
- Adopted in `families.yaml` for all three (mrpB, mrpF, mrpG); see each
  family's own entry for the full writeup. mrpG's redesign wasn't backed
  by mrpB/mrpF's volume of live-UniProt-sampling evidence going in (a
  reasonable bet based on its identical `pfam_model` +
  `pfam_ga_review_needed` shape, not independently confirmed the same
  way) -- but the 7-family benchmark result validates the outcome
  directly regardless.
- gshB's negative pool not touched -- confirmed structurally different
  (see above), the mrp-style redesign doesn't apply. Its `pfam_model` was
  a separate, real problem though -- see below.

### gshB: pfam_model removal recovers the biggest single-family F1 win on this branch

Panel re-ranked after the murB/mrpB/mrpE/mrpF/mrpG fixes above to find the
next priority. gshB stood out with a distinctive profile unlike any of the
mrp genes: HMM precision was already excellent (0.968) -- confirming its
negative pool (63% purity-flagged, the worst in the panel by that raw
number) was never actually the problem, consistent with the real minority
diversity found in it above -- but HMM recall was catastrophic (0.307) on
a large positive pool (7,756), ruling out data-thinness as an explanation
too.

- **Root cause**: gshB uses `pfam_model: PF02951`, the same setup already
  found to hurt mrpE's recall. gshB is a ~310aa protein; Pfam's GA cutoff
  is calibrated for full-length curated-sequence matches, not the short
  fragments 6-frame-translated simulated reads produce, so real positive
  reads were legitimately scoring below it. The reasoning that originally
  justified keeping `pfam_model` here ("a contaminated pool couldn't
  produce a reliable local calibration anyway, so deferring to Pfam's GA
  cutoff is the reliable path") is the identical logic already disproven
  on mrpE -- it validated the calibration-status flag, not the real F1.
- **Tested as a 4-family subset** (gshA, gshB, gshF, gor via
  `make_family_subset.py`), only gshB's `pfam_model` removed, siblings
  left unchanged as a cross-family check -- none share PF02951 with gshB,
  so little competition risk was expected going in, and none appeared:
  gshA/gshF/gor's numbers moved by single digits (noise) between the two
  runs, both before and after. DIAMOND unaffected in all four families
  (`pfam_model` only touches the HMM path).
- **Result**: HMM recall 0.307->0.827, precision 0.968->0.942 (small,
  acceptable cost), F1 0.466->0.880 -- the largest single-family F1 gain
  found on this branch, bigger than mrpE's own version of the same fix
  (0.19->0.704, which combined a negative-pool redesign with the
  `pfam_model` removal; gshB got most of its win from `pfam_model` alone).
  Adopted in `families.yaml`; negative_query and negative_pfam unchanged.

### opuBA: paralog best-hit competition root-caused; RefSeq expansion tried and reverted

Panel re-ranked again after gshB. opuBA stood out with yet another distinct
profile: HMM/DIAMOND precision already perfect (1.0/1.0), recall
catastrophic (0.099/0.29) -- ruling out both the mrp-style negative-pool
mechanism and the gshB-style `pfam_model` mechanism (opuBA uses neither).

- **Root cause, traced directly through `gene_counts.tsv`/
  `hmmscan.tblout`** rather than assumed: `11_compute_metrics.py`'s
  `best_hit_per_read()` assigns each read to whichever family scores
  highest across the *entire* combined database, with no cutoff involved
  in that count. Checked where opuBA's own true-positive reads actually
  go: the majority win best-hit against **opuCA's** model instead of
  opuBA's own (HMM: 140/230 hit reads to opuCA, 58 to opuAA, only 32 to
  opuBA; DIAMOND shows the same skew). Same category of problem as
  opuBC/opuAC (close-paralog competition) -- but unlike opuBC, no
  `decoy_from_negatives` setting caused it; opuAA/opuBA/opuCA are
  ordinary, always-co-present families, so there's no config to revert.
  opuBA's own UniProt-curated pool (157) is roughly half opuAA's (325)
  and opuCA's (358), a plausible reason it loses the competition more.
- **RefSeq expansion attempted** (same mechanism as 9 already-enabled
  families), verified with the established discipline first: a 500-entry
  live-NCBI sample showed 56.9% novelty at 90% identity (healthy) and a
  real but minor collision (4.8% explicitly titled "...ProV", a
  different gene) -- small enough to handle with a new, narrowly-scoped
  `refseq_exclude_products` filter added to `01e_add_refseq_positives.py`
  (substring-matches a RefSeq entry's product description, skips it at
  fetch time) rather than excluding the family outright the way
  mrpC/mrpF were for their own, larger collisions.
- **Verification sample turned out to be built from the wrong query,
  though** -- caught only after actually running it: `01e`'s own
  `build_query()` restricts to `srcdb_refseq[Properties]` (NCBI's curated
  non-redundant subset), matching only 236 entries for opuBA, not the
  ~2,500 the broader verification query found. Real fetch: 233 sequences
  appended, but CD-HIT clustering (90% identity, the same step every
  family's pool goes through) collapsed 232 of them into existing
  clusters -- only 1 was genuinely novel. The curated RefSeq subset
  already overlaps almost entirely with UniProt's own opuBA population;
  the broader verification population looked more novel specifically
  because it included redundant/less-curated entries CD-HIT discards
  regardless.
- **Result** (3-family opuAA/opuBA/opuCA subset build+benchmark,
  `make_family_subset.py`): DIAMOND flat (F1 0.459->0.463, one extra
  training sequence), **HMM regressed** (F1 0.289->0.167) -- consistent
  with a single added sequence perturbing the alignment/model build as
  noise, not supplying real signal. **Reverted** -- opuBA's `families.yaml`
  entry back to its original queries, with the full investigation
  documented there. The `refseq_exclude_products` mechanism itself is
  kept in `01e_add_refseq_positives.py` (generic, tested, currently
  unused by any family) for the next case that needs it, same reasoning
  the opuBC decoy-dedup fix was kept in `01_fetch_refs.py` after its own
  revert.
- **Lesson for next time**: verify RefSeq novelty/collision against
  `01e`'s own `build_query()` output directly, not an approximation of
  it -- the `srcdb_refseq[Properties]` restriction materially changes
  both population size and its overlap with UniProt.
- opuBA's real problem (paralog best-hit competition) remains unsolved,
  same open-question status as opuBC/opuAC.

### cspA: last remaining pfam_model family, small confirmed win

Audited the panel for any other family still using `pfam_model` after the
mrpE/gshB fixes above -- only cspA. Checked rather than assumed it would
behave the same way: precision already perfect (1.0/1.0 both methods),
recall showing the same shape of gap (HMM 0.627 vs. DIAMOND 0.681) as
mrpE/gshB had, and `qc_scorecard.tsv` showing the cleanest possible
starting point for this fix (`n_negative_scored=0`, zero overlap between
cspA's negative pool and its positives -- no redesign needed, just the
cutoff source).

- Isolated single-family build+benchmark (`make_family_subset.py`,
  `pfam_model` removed, `negative_query` unchanged) confirmed a real but
  modest gain: HMM F1 0.771->0.783 (recall 0.627->0.643), precision
  unchanged. Adopted -- no downside found.
- Notably smaller than mrpE's or especially gshB's version of the same
  fix. Plausible reason: cspA's positive set is deliberately broad by
  design (captures every csp paralog across all of Bacteria, not a single
  ortholog group, per the family's own Task 1 spec), so the custom-built
  HMM ends up close to Pfam's own long-established, well-curated CSD
  model anyway -- there wasn't much slack between the two cutoff sources
  to recover, unlike mrpE/gshB's single-ortholog-group models.
- `pfam_model` is now unset for all 43 families -- this lever is fully
  exploited across the panel, nothing left to audit here.

### mrpA/mrpD: HMM underperformance root-caused, a genuine mixed problem (not yet fixed)

Investigated why HMM trails DIAMOND by a wide margin for both families
(mrpA: HMM F1 0.382 vs. DIAMOND 0.736; mrpD: HMM F1 0.428 vs. DIAMOND
0.64) -- distinct from every fix above, confirmed by checking rather than
assuming it was the same mechanism.

- Best-hit competition is NOT the driver here (unlike opuBA): mrpA/mrpD
  correctly win ~94% of their own true-positive reads' best-hit calls.
- Traced the false positives directly (cross-referencing every other
  family's own `hmmscan.tblout` for hits landing on mrpA's/mrpD's model):
  they come overwhelmingly from **each other's own negative sets**
  (mrpA's FPs: 5,078 from mrpD.negative, 4,780 from mrpB.negative;
  mrpD's FPs: 6,593 from mrpA.negative, 1,965 from mrpB.negative), not
  from unrelated families scattered across the panel. The mrpB
  contribution is likely stale -- that `v9` data predates this branch's
  own mrpB fix (mrpB's old negative pool was ~68% genuine mrpA+mrpB
  fusion sequences, which would legitimately score well against mrpA's
  model; not re-verified with a fresh build).
- A 200-entry live UniProt sample of the shared PF00361 negative pool
  (mrpA/mrpD's own domain) shows this is a genuine **mix**, unlike every
  fix above: ~56% are real, different, legitimately-hard-negative
  proteins (NADH-quinone oxidoreductase subunits L/M/N from Complex I,
  and DabB, both structurally related via a shared broad fold but
  functionally distinct -- exactly the "generic Complex I NuoL/M/N-type
  subunits" caution already written into mrpA's own family description),
  ~15% are genuinely mislabeled mrpA/mrpD orthologs under bare locus tags
  (the same structural-contamination pattern already fixed elsewhere),
  and the rest ambiguous from title alone.
- **Not fixed yet, deliberately**: a full random-background swap (the
  mrpB/E/F/G playbook) would work for the mislabeled-positive ~15% but
  would also discard the legitimate ~56% Nuo/DabB hard-negative content
  that mrpB/E/F/G's negative pools never had -- real risk of trading a
  precision win for a specificity loss against genuinely confusable
  Complex I subunits in real samples, unlike every prior fix on this
  branch. A more surgical fix (exclude only the specific
  mislabeled-positive product-description pattern from each other's
  negative pools, keep the Nuo/DabB content) was flagged as the more
  promising untried idea -- tried next (below), and it didn't work.

### mrpA/mrpD: surgical exclusion tried, made both methods worse -- reverted

Followed up on the surgical-exclusion idea flagged above, using the same
methodology just proven out on mscS.

- Built and live-verified an exclusion query targeting only the confirmed
  mislabeled-positive pattern (`NOT (protein_name:antiporter AND
  (protein_name:"subunit A" OR protein_name:"subunit D" OR
  protein_name:MnhA OR protein_name:MnhD))`) -- a 200-entry check of what
  it actually matched came back almost entirely genuine mislabeled mrpA/
  mrpD content (various subunit A/D and MnhA/MnhD phrasings, plus CPA3
  family classification), no legitimate Nuo-only entries caught by
  mistake. ~18.5% of the shared PF00361 pool excluded (13,174/71,137).
- Tested directly on mrpA's and mrpD's real `negative_query` (not a
  calibration-only side channel like mscS's `hmm_negative_query` -- the
  false positives here were confirmed to be real best-hit-competition
  losses, so a real fix needed to move the actual reads/model, and
  `hmm_negative_query` wouldn't have been the right tool for that).
  7-family mrpA-G subset build (`make_family_subset.py`), baseline vs.
  surgical.
- **Result: worse across the board**, both families, both methods --
  mrpA DIAMOND F1 0.695->0.666 (fp 1048->1228), mrpA HMM F1 0.383->0.359
  (fp 2980->3361), mrpD DIAMOND F1 0.629->0.601 (fp 1739->1924), mrpD HMM
  F1 0.419->0.386 (fp 3364->3657). False positives went UP, not down, in
  every row.
- **Why, most likely**: unlike mrpB/mrpE/mrpF/mrpG's full random-
  background swap -- where ~95-99% of the pool was contamination, so
  replacing all of it gave calibration an unambiguous, clean gap to
  anchor on -- this pool was already a genuine mix with substantial
  (~67%) legitimate hard-negative content. Removing only the highest-
  scoring, most positive-like subset (the mislabeled ~18.5%) while
  leaving the more distant Nuo/DabB majority in place plausibly let
  DIAMOND/HMM's score-distribution-based calibration drift toward a more
  permissive threshold fit to the now-easier remaining pool -- the same
  calibration-drift mechanism already seen on mrpB's fusion-fix history
  and mscS's own surgical test, but here hurting both methods at once
  instead of trading one for the other.
- **Not adopted.** Documented directly in mrpA's `families.yaml` entry so
  this isn't retried blindly. mrpA/mrpD's HMM underperformance remains
  genuinely open -- the lesson generalizes: partial/surgical cleaning of
  an already-mixed negative pool is not a smaller-scale version of the
  full-replacement fix that worked on near-totally-contaminated pools;
  it's a different, apparently unfavorable trade.

### mscS: surgical negative-pool exclusion, and the per-method-cutoff mechanism it uncovered

Investigated the same way as mrpA/mrpD above -- traced false positives
(negligible cross-panel contamination, confirming they come from within
mscS's own negative pool) and sampled it live: a genuine mix of real named
paralogs (kefA, ybdG, ykuT) and a large bare-locus-tag population whose
product description is literally "Mechanosensitive ion channel protein
MscS" -- mislabeled positives.

- **Surgical exclusion query** built and tested: `xref:pfam-PF00924 NOT
  gene:mscS* AND (protein_name:KefA OR protein_name:YbdG OR
  protein_name:YkuT OR NOT protein_name:MscS)` -- keeps the named
  paralogs, drops the bare-tag entries (54,533 -> 2,757 candidates).
- **Full swap first** (both DIAMOND and HMM calibrated + benchmarked on
  the surgical pool): HMM improved (F1 0.481->0.546) but **DIAMOND
  regressed** (F1 0.546->0.5) -- the same calibration-drift mechanism
  already documented on mrpB's fusion-fix history recurring in a new
  place: `08b_calibrate_diamond_cutoffs.py` sets its cutoff from the
  negative set's own distribution, and removing the (conceptually wrong
  but score-similar) bare-tag entries loosened it.
- **New mechanism built to avoid that trade**: `hmm_negative_query`, an
  opt-in families.yaml field (implemented in `01_fetch_refs.py` +
  `05_build_hmms.sh`) that lets a family's HMM calibration draw from a
  different negative pool than DIAMOND's -- the negative-side mirror of
  `08_build_diamond_db.sh`'s existing DIAMOND-specific positive-side split
  (`_study`/`_refseq` exclusion). Set on mscS to the surgical query above.
- **Verified the deployed mechanism directly, not assumed** -- and found
  it doesn't do what was expected: re-tested after wiring it in, result
  was byte-identical to the unfixed baseline. Root cause: `11_compute_
  metrics.py`'s summary.tsv (the metric every F1 number in this file
  comes from) is pure best-hit-per-read competition on raw hmmscan bit
  scores and never applies the calibrated GA cutoff at all (confirmed by
  reading `load_hmm_calls`/`best_hit_per_read`/`label_reads` directly).
  `hmm_negative_query` only feeds the GA-cutoff calibration, so it cannot
  move this metric for a family whose model and benchmark reads are
  otherwise unchanged.
- **Checked what production actually consumes** before deciding whether
  to keep it: `build_profile_cascade_config()`'s own docstring confirms
  osmotool's `profile`-mode DIAMOND+HMM cascade uses HMM's short-read-
  swept `best_threshold_summary.tsv` value, explicitly NOT the GA cutoff
  ("these are short reads, not full-length ORFs") -- so `profile` mode is
  provably unaffected by this fix either way. `annotate` mode, by
  contrast, does use the GA cutoff directly, on full-length ORFs -- the
  same data type this calibration actually scores. There, the fix is
  real: calibration-time whole-protein separation (`f1_at_cutoff` in
  `cutoff_manifest.tsv`) went from 0.697 to 0.815, cutoff shifted
  221.9->214.1 bits.
- **Decision: kept**, for `annotate` mode specifically -- a correctly-
  targeted, low-risk data-quality fix (removes genuine mislabeling from
  what HMM calibrates against) with no effect on `profile` mode either
  way. Not independently benchmark-verified for `annotate` mode's real
  accuracy, since this pipeline has no whole-ORF-level benchmark the way
  it has one for short reads -- flagged honestly in mscS's own
  `families.yaml` entry rather than overstating it as a measured win.
- **Reusable takeaway for the rest of this panel**: `hmm_negative_query`
  now exists as a real, tested mechanism for the next family that needs
  DIAMOND and HMM to disagree on what a hard negative is -- but it should
  only be reached for with a clear-eyed understanding that it changes
  `annotate`-mode behavior, not this repo's own short-read benchmark
  numbers or `profile`-mode production behavior.

### opuBB: fourth confirmed instance of the opuBC/opuAC/opuBA close-paralog problem

Final panel-audit check before considering the branch release-ready --
opuBB was the one remaining family with an unexamined "good precision,
poor recall" profile (DIAMOND 0.921/0.257) matching opuBA's own signature
before its root cause was found.

- Traced directly (`gene_counts.tsv`/`hmmscan.tblout`, not assumed): the
  majority of opuBB's own true-positive reads win best-hit against
  **opuCB's** model instead of its own (DIAMOND 62 vs. 35.5; HMM 119 vs.
  71) -- the identical close-paralog best-hit-competition mechanism
  already root-caused on opuBC/opuAC/opuBA, just in the permease-subunit
  trio (opuAB/opuBB/opuCB, PF00528) instead of the ATPase one
  (opuAA/opuBA/opuCA, PF00571).
- Same population pattern too: opuBB's UniProt-curated pool (187) is the
  thinnest of its trio (opuAB 342, opuCB 256) -- opuAB itself performs
  fine (F1 0.82-0.88), consistent with thin-pool competitive loss being
  the driver rather than a domain-wide problem.
- RefSeq expansion not attempted -- checked the curated population size
  first given opuBA's own attempt already failed from near-total
  redundancy with UniProt after clustering: opuBB's curated RefSeq
  population (141) is smaller than opuBA's was (236), making the same
  outcome likely enough not to spend a full build+benchmark test
  confirming it.
- Documented in opuBB's own `families.yaml` entry. Same open-question
  status as opuBC/opuAC/opuBA -- no fix available with tools tried so
  far. With this check done, every family flagged during this session's
  full-panel re-ranking has now been either fixed or root-caused and
  documented as a genuine dead end.

### Full 43-family `v8` rebuild -- panel-wide validation of every fix on this branch

Built and benchmarked the full panel (`releases/v8`) to confirm the 8
fixes above hold up at real full-panel scale, not just in the isolated/
small-group subset builds each was validated in individually -- the same
"validate the whole branch before shipping" discipline this project used
for `v7`/`v9`. Named `v8` since the previous local `v8` build (from an
earlier session, `cluster-negatives-cdhit` branch) was never archived/
published (no `v8.tar.gz`, unlike `v5`/`v7`) -- moved aside to
`releases/v8-prerelease-backup` rather than deleted, in case it's still
wanted.

- **All 8 fixed families held up cleanly**, several with dramatic gains
  that match (and in some cases modestly exceed) their isolated-test
  numbers: mrpB HMM F1 0.249->0.779 (fp 713->0), mrpE HMM F1
  0.196->0.717 (fp 2311->0), mrpF HMM F1 0.271->0.718 (fp 790->0), mrpG
  HMM F1 0.484->0.619 (fp 758->0), gshB HMM F1 0.467->0.878 (fp mostly
  unchanged, precision held at 0.941), murB HMM F1 0.705->0.738, cspA
  HMM F1 0.756->0.789, mscS DIAMOND F1 0.53->0.541 / HMM F1
  0.465->0.501. Precision stayed at or near 1.0 for the mrp
  genes/murB/cspA even with all 42 other families' sequences competing
  in the same DIAMOND/HMM database -- confirms the isolated/small-group
  validation wasn't hiding a cross-family regression the way mrpB's own
  fusion-fix history warned it could.
- **Documented dead ends re-confirmed stable, not a testing artifact**:
  mrpA/mrpD, betB, and opuBC/opuAC/opuBA/opuBB all moved by only normal
  full-panel-composition drift (consistent with this project's own
  documented noise level for untouched families), nothing alarming --
  their open-question status stands.
- **Panel-wide (read-volume-weighted)**: DIAMOND F1 0.827->0.832, HMM F1
  0.693->0.722. A modest aggregate shift is the correct outcome given
  only 8 of 43 families changed and the other 35's much larger combined
  read volume dilutes the average -- the real story is in the 8
  families' own numbers above, not the panel-wide one.
- `releases/v8` is a local build only (`releases/` is gitignored) -- no
  further action needed to preserve it beyond this changelog entry,
  same as every prior release in this project's history.

### Status at time of writing

- RefSeq merge (`01e`) and the per-method reference split
  (`08_build_diamond_db.sh`) are both implemented and validated at three
  levels: sequence-merge correctness, single-family build+benchmark (4
  families), and now a full 43-family panel rebuild (`v9`, above) --
  confirming the split's benefit holds at the scale this will actually
  ship at, not just in isolation.
- UniProt/RefSeq provenance tracking fixed, Docker image verified working
  and its osmotool pin brought current, `.dockerignore` added, and
  murB/mrpE/mrpB/mrpF/mrpG's negative-pool fixes plus gshB's `pfam_model`
  removal all validated (mrpE alone in isolation then re-validated with
  mrpB/F/G as a 7-family group; the other three validated directly as a
  group; gshB validated as its own 4-family gsh group) and made permanent
  in `families.yaml` -- all committed on this branch. opuBC/opuAC
  re-examined and confirmed to be a different, already-exhausted problem
  (close-paralog competition, not negative-pool mislabeling) -- not
  touched. opuBA root-caused to the same close-paralog-competition
  category (loses best-hit against opuCA); a RefSeq expansion attempt was
  tried and reverted after the real fetch turned out almost entirely
  redundant with UniProt (only the verification sample looked novel, from
  a query that didn't match production). opuBA's underlying problem
  remains open, same status as opuBC/opuAC. cspA's `pfam_model` also
  removed (small, real HMM recall gain, no downside) -- no family in the
  panel uses `pfam_model` anymore. mrpA/mrpD root-caused to a genuine
  mixed problem (real Complex-I-subunit hard negatives alongside
  mislabeled-positive contamination in their shared PF00361 pool); a
  surgical exclusion of the mislabeled subset was tried and made both
  methods worse for both families (calibration drift, same mechanism as
  mrpB's fusion-fix history) -- reverted, documented, remains genuinely
  unfixed with what's been tried. betB confirmed unfixable with current tools (promiscuous domain, no
  narrower anchor, no single identifiable contamination pattern). mscS
  fixed for `annotate` mode via a new `hmm_negative_query` mechanism
  (per-family HMM-only negative pool, mirroring `08_build_diamond_db.sh`'s
  DIAMOND-specific positive-side split on the negative side) -- verified
  it does NOT affect `profile` mode or this repo's own short-read
  benchmark numbers (both use best-hit/swept-threshold logic that bypasses
  the GA cutoff this field feeds), kept anyway for its real, if
  unbenchmarked, `annotate`-mode calibration improvement. opuBB confirmed
  as a fourth instance of the opuBC/opuAC/opuBA close-paralog problem
  (loses best-hit against opuCB) -- same open-question status, documented,
  not touched. This closes out the full-panel re-ranking done this
  session: every family flagged as underperforming has now been either
  fixed or root-caused and confirmed as a genuine, documented dead end.
  Full 43-family `v8` rebuild confirms all 8 fixes hold up at real
  panel scale (see that section above) -- this branch is release-ready.
- Committed and pushed to `add-refseq-positives` (PR #15, still open).
- `v9` (and the murB/mrpG/otsA/mscS/otsB single-family test releases) are
  local builds only (`releases/` is gitignored) -- no further action
  needed to preserve them beyond this changelog entry.
