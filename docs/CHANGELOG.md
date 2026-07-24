# Changelog

Chronological summary of osmo_refdb's development, from the initial
pipeline build through the current negative-pool quality follow-up work
on branch `fix-numbered-paralog-gaps`. Full rationale and evidence for
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
