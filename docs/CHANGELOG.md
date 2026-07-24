# Changelog

Summary of the negative-pool quality follow-up work on branch
`fix-numbered-paralog-gaps` (built on top of PR #6's Task 1 gene panel and
PR #7's mrpF/mrpG PhaF/yufB synonym fix). Full rationale and evidence for
each item lives in `families.yaml`'s per-family `description` fields and
`README.md`; this file is a chronological summary, not the source of truth.

## Investigation (Context)

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

## Phase 1 — Cleanup

Discarded the `max_negative_override` experiment's code changes (kept the
finding, not the code). No behavior change.

## Phase 2 — Numbered-paralog gap fix (issue #8, PR #9)

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

## Phase 3 — Structural negative-pool contamination (documentation)

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

## Phase 4 — proX/opuAC/opuBC/opuCC decoy conversion (tried, disproven, reverted)

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

## Status at time of writing

- **Phase 1, 2**: committed on `fix-numbered-paralog-gaps`; issue #8 filed;
  PR #9 open against `fix-mrpF-mrpG-synonym-gap`.
- **Phase 3**: implemented and validated (`yaml.safe_load`, 43 families,
  no dupes) in an isolated git worktree (`.worktrees/phase3`, branch
  `phase3-structural-docs`) to avoid editing `families.yaml` while the
  `v7` full-rebuild Docker container has it mounted live. Pending `v7`'s
  completion before merging back and committing alongside Phase 2's sweep.
- **Phase 4**: disproven and reverted (see above) — no changes pending
  merge for this phase; the worktree's `families.yaml` for
  proX/opuAC/opuBC/opuCC now matches PR #6's original design again.
- **Phase 5**: reframed by Phase 4's finding — trkH/ktrB/ktrD share a
  Pfam accession the same way proX/opuAC/opuBC/opuCC did, so the decoy
  option there is now disfavored by default rather than a neutral
  re-evaluation. Not yet started.
- **Phase 6** (qc_scorecard-as-gate documentation for `pfam_model`): not
  yet started.
