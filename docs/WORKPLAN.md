# Workplan: osmo_refdb negative-pool quality follow-ups

Saved copy of the active plan, in case the working session crashes before
this lands as commits. See `docs/CHANGELOG.md` for a narrative summary of
what's already done; this file is the working plan itself.

**Status**: Phases 1-4 are done and committed on `fix-numbered-paralog-gaps`
(commits `69f80a4`/`6d03487`/`da5caf2`/`837e600`, pushed, PR #9 and issue
#8 updated). Phase 5 is done and committed on `phase5-trk-ktr-scope`
(commit `8b49d2d`, pushed, PR #11 open against `fix-numbered-paralog-gaps`).
Only **Phase 6** is still open — see that section below. The
crash-recovery notes that used to live here (git worktree state, pending
`v7` build, uncommitted changes) are no longer relevant now that
everything has landed; kept only in `git log` if needed for reference.

## Resuming after a crash — check these first

1. **Git state**: branch `fix-numbered-paralog-gaps` should have two
   commits on top of `fix-mrpF-mrpG-synonym-gap`: `69f80a4` (original
   6-family fix) and `6d03487` (full 43-family sweep + Phase 3 docs +
   Phase 4 finding). If `6d03487` is missing, Phases 2-4 need redoing —
   full content is described below and was also in `git log`'s commit
   message if the commit exists but wasn't pushed.
2. **GitHub**: issue #8 (numbered-paralog gap bug class) and PR #9
   (`fix-numbered-paralog-gaps` → `fix-mrpF-mrpG-synonym-gap`) are open,
   both updated with the final findings. If a crash happened before the
   `gh pr comment`/`gh issue comment` calls landed, the PR/issue may still
   only show the original 6-family framing even though the commit itself
   is complete and pushed — just re-run the comment calls (text is in
   `git log`'s history of this conversation, or can be re-derived from
   `docs/CHANGELOG.md`).

## Context

PR #6 (Task 1 gene panel) and PR #7 (mrpF/mrpG PhaF/yufB synonym fix)
surfaced a family of related data-quality issues in how per-family
hard-negative pools are built. Evidence was pulled from the last clean
`v6` build (`releases/v6/refs/negative_purity_manifest.tsv` +
flagged-sequence gene lookups via UniProt) and one experiment was run
live before writing this plan — its result **disproves** the most obvious
fix and reshapes the plan materially:

- Raising mrpF/mrpG's negative cap 1000→4000 (`max_negative_override`)
  gave 4-5x more *clean* negatives after purity filtering (413→2120/1621)
  but made benchmark F1 **worse**, not better (mrpF DIAMOND F1
  0.329→0.109, mrpG 0.535→0.319). More volume from the same narrow Pfam
  pool dilutes with borderline/confusable cases rather than adding real
  discriminative signal. Reverted, not merged — finding documented in
  Phase 3 instead.
- Checking flagged sequences for mrpB (23%), mrpE (21%), mrpC (26%), gshB
  (63%), otsA (58%), mazG (64%), mscL (72%), and most of trkH
  (27%)/ktrB (14%) shows the **same root cause** every time: the
  "negative" pool drawn from the same narrow Pfam family is mostly real
  orthologs of the target gene that UniProt never assigned a curated gene
  symbol to (bare locus tags) — not a missing-synonym bug, not a
  threshold-calibration artifact. The flat 70% purity threshold is not
  the problem; if anything the `max_negative_override` result suggests it
  may already be too permissive, not too strict.
- Two genuinely different, cheap-to-fix gaps turned up along the way,
  same *class* of bug as PhaF/yufB (documented-alias-style) but a new
  pattern (numbered paralogs): missing `<gene>1`/`<gene>2` UniProt gene
  symbols — a genuine undercount of true positives, not a
  negative-contamination problem. See Phase 2.
- proX/opuAC/opuBC/opuCC show *low* purity-flag rates (the cross-exclusion
  added in PR #6 works) but still poor precision — genuine cross-paralog
  confusion between four real, distinct genes, the exact shape
  `decoy_from_negatives` (proven on betL) exists to fix. See Phase 4.

## Phase 1 — Clean up the disproven experiment — DONE

Discarded the uncommitted `max_negative_override` edits on what's now
`fix-numbered-paralog-gaps`. Finding folded into Phase 3's documentation
instead of the code.

## Phase 2 — Numbered-paralog gap fix — DONE (committed, PR #9)

Branch `fix-numbered-paralog-gaps` off `fix-mrpF-mrpG-synonym-gap`.

- Started as a 6-family fix (murB, trkA, trkH, ktrA, ktrB, ktrD), each
  numbered variant verified per-accession against UniProt to confirm the
  right sibling family.
- Subset-tested against `v6`: **not** a uniform win. murB/trkH improved;
  ktrA/ktrB/ktrD's swings were consistent with noise from small held-out
  test sets (ktrD's DIAMOND cutoff calibrated on only 9 test positives);
  **trkA regressed** (DIAMOND F1 0.920→0.879).
- trkA's regression was root-caused, not just reverted: of 372 UniProt
  `trkA1`/`trkA2` entries, 82% (median 223aa) are short single-domain
  fragments already excluded by the existing length-outlier filter; the
  68 that survive are real, taxonomically diverse full-length orthologs
  whose diversity widens the DIAMOND reference cloud and costs
  specificity. Kept anyway (real genes, not an artifact) and documented
  in `families.yaml`. Same underlying principle as the
  `max_negative_override` disproof — broader coverage isn't free, it
  trades against discriminative sharpness, just triggered from the
  positive side here instead of the negative side.
- Filed issue #8, opened PR #9 (`fix-numbered-paralog-gaps` →
  `fix-mrpF-mrpG-synonym-gap`). This part is **committed**
  (`69f80a4 Add numbered-paralog gene symbols to
  murB/trkA/trkH/ktrA/ktrB/ktrD queries`).
- **Expanded to a full systematic sweep** of all 43 families after
  discovering this bug class is near-universal: every single-token gene
  symbol in the panel has real, curated numbered variants in UniProt,
  from a handful (mazG: 9) to hundreds (nhaA: 304/389, galE: 678/530).
  Verified as real curated data (not a search-engine artifact) via direct
  UniProt lookups (e.g. `A8LVS8 nhaA1` in *Salinispora arenicola*). `cspA`
  is structurally exempt (anchored on the bare Pfam accession, not a gene
  symbol, by design).
  - 16 "singles" families (ectA/B/C, betL cluster, kdpA, nhaA, proP,
    otsA/B, mscL/S, galE, mazG, betA/B, gor, gshB): mechanical
    `<sym>1`/`<sym>2` addition.
  - 6 shared-Pfam clusters requiring per-accession cross-family
    verification (to rule out cross-contamination): proX/opuAC/opuBC/opuCC
    (PF04069), opuAA/opuBA/opuCA (PF00571), opuAB/opuBB/opuCB (PF00528),
    mrpA/mrpD (PF00361), mrpB/mrpC/mrpE/mrpF/mrpG (their own individual
    Pfam accessions, mnh/PhaF/yufB aliases numbered too), gshA/gshF
    (PF04262). All checked clean.
  - Given 36 of 43 families changed, treated as a full rebuild
    (`releases/v7`) rather than a subset test, per the project convention
    of using a new version number for any full rebuild.
  - **Verified against `v6`**: 22-25 families improved per method vs.
    17-18 regressed above noise level (net positive per-family), but
    read-volume-weighted precision/recall is essentially flat (DIAMOND F1
    0.823→0.821, HMM 0.684→0.684) — a completeness fix, not a performance
    win. Full numbers in `docs/CHANGELOG.md`. All regressions >0.03 F1
    documented in their family's `families.yaml` entry.
  - **Committed** as `6d03487` on `fix-numbered-paralog-gaps`, pushed,
    issue #8/PR #9 updated with the final findings.

## Phase 3 — Document the structural finding — DONE, committed

Committed together with Phase 2's sweep in `6d03487` (was implemented in
a git worktree during development, since removed — no longer relevant).

- Updated `families.yaml` entries for mrpB, mrpC, mrpE, gshA, gshB, gshF,
  otsA, mazG, mscL: purity contamination here is structural (Pfam family
  mostly coextensive with the target gene; UniProt's automatic/unreviewed
  annotations don't assign gene symbols to most members), confirmed via
  direct sequence lookup, and not fixable by raising
  `max_negative_override` (disproven on mrpF/mrpG) or widening
  `positive_query` (disproven on trkA).
- **`pfam_model` vindication**: for mrpB, mrpE, gshB (all flagged
  `pfam_ga_review_needed` in `cutoff_manifest.tsv`), this structural
  finding validates the original `pfam_model` adoption for a second,
  independent reason — a locally-calibrated negative set couldn't have
  been trusted here anyway.
- **Production disposition, checked individually, not assumed uniform**:
  of the four families with no `pfam_model` fallback (otsA, mazG, mscL,
  gshA), `profile_cascade.tsv`'s DIAMOND→HMM fallback doesn't rescue any
  of them (HMM calibration draws on the same contaminated pool — all four
  carry `hmm_status=overlapping_distributions_f1_calibrated`, not an
  independently-calibrated status). Checked DIAMOND precision per family
  rather than assuming they all need the same treatment:
  - otsA (0.543) and mrpC (0.479, found during this check to be worse
    than otsA despite not being in the original four-family list) — set
    `scope: annotate_only`.
  - mazG (0.870) — borderline, kept in profile mode; already carries its
    own biological-confidence caveat in `families.yaml`.
  - mscL (0.927), gshA (0.905) — fine in practice, calibration still
    found a workable cutoff despite the noisy pool. No scope change.
- Added a new README section, "Structural negative-pool contamination:
  when more data can't fix it," proofread for line-wrap issues (per user
  instruction to always check README diffs for this).

## Phase 4 — proX/opuAC/opuBC/opuCC decoy conversion — TRIED, DISPROVEN, REVERTED

Implemented and subset-tested (`releases/v-test-phase4`, 4-family build)
in `.worktrees/phase3`; test artifacts deleted after use, nothing to
resume. **Do not redo this without a specific reason to expect a
different outcome** — see the finding below.

- Redesigned all four families' `negative_query` to target the other
  three specifically as the decoy source (mirroring betL vs. betT/caiT)
  instead of drawing from the broad, anonymous PF04069 pool, with
  `decoy_from_negatives: true` set on all four. Added a dedup check to
  `pipeline/08a_build_decoy_refs.py` first (`load_all_positive_accessions()`,
  syntax-checked, **kept in the codebase** — reads every family's
  `refs/*.positive.train.faa` accessions into one exclusion set and
  filters decoy candidates against it before writing, so a decoy set
  can't contain a sequence that's also a real positive reference
  elsewhere; still a valid general safeguard for betL and any future
  `decoy_from_negatives` use).
- **Subset-tested before merging — result was a catastrophic DIAMOND
  recall collapse for all four families**, not the modest recall cost
  anticipated: opuBC recall dropped to exactly 0.0 (of 240 test reads,
  every one that aligned was attributed to a sibling's decoy label
  instead of opuBC's own reference); proX/opuAC/opuCC all fell under
  0.08. Root-caused via `opuBC.positive.gene_counts.tsv`.
- **Mechanism**: unlike betT/caiT (only superficially score-confusable
  with betL, genuinely distinguishable at the sequence level),
  proX/opuAC/opuBC/opuCC are close enough paralogs of *each other* that
  adding them back into the searchable pool as decoys lets them
  systematically outcompete a family's own thin reference set (opuBC in
  particular has only 93 UniProt members to begin with). The dedup fix
  addressed the specific collision it was designed for (same accession as
  both decoy and positive reference) but not this failure mode
  (genuinely similar, non-identical sequences winning on merit) — a
  useful distinction for judging any future `decoy_from_negatives`
  proposal.
- **Reverted**: proX/opuAC/opuBC/opuCC back to PR #6's cross-exclusion-only
  design, committed as part of `6d03487` alongside Phases 2/3 (the
  finding documented in each family's description, not a functional
  change — these four families' queries match what was already in
  production). Their remaining precision problem (genuine cross-paralog
  confusion) is an open question without a known fix, not a pending
  decision.

## Phase 5 — Re-evaluate trkH/ktrB/ktrD decoy need — DONE (PR #11)

Checked `v7`'s real numbers rather than assuming either fix helps:

- Numbered-paralog fix (Phase 2) barely moved purity-flag rates (trkH
  27%→29%, ktrB 14.5%→15.5%, ktrD 2.4%→2.8% flagged, v6 vs v7) — same
  structural pattern as Phase 3's families, not a missing-symbol gap.
- `decoy_from_negatives` ruled out per Phase 4's finding — trkH/ktrB/ktrD
  share PF02386 the same close-paralog way proX/opuAC/opuBC/opuCC shared
  PF04069, the exact clique shape that catastrophically backfired.
- No `pfam_model` fallback possible (removing it fixed the byte-identical-
  HMM bug on this same trio).

Both DIAMOND (0.265–0.293) and HMM (0.121–0.288) precision are low with
no independent cascade rescue. Set `scope: annotate_only` on all three,
same criteria as otsA/mrpC — confirmed explicitly with the user first
since this removes the entire Trk/Ktr K+-uptake system from
`osmotool profile`'s output, a bigger call than extending scope on 1-2
families. Committed `8b49d2d` on branch `phase5-trk-ktr-scope`, pushed,
PR #11 open against `fix-numbered-paralog-gaps`. Also backfilled mrpC's
missing README documentation from the prior commit.

## Phase 6 — qc_scorecard as a mandatory gate for pfam_model — NOT STARTED

New branch `pfam-model-gate-docs`, lowest urgency, independent of the
above.

- Update `README.md`'s `pfam_model` section: state plainly that
  `domain_architectures` count is a screening heuristic, not sufficient
  on its own (5 of 7 adoptions this round got flagged
  `pfam_ga_review_needed`) — `qc_scorecard.tsv` review after the first
  build is a required step before trusting a `pfam_model` adoption, not a
  follow-up nicety.
- Consider (lower priority, only if time allows): a small enhancement to
  `06b_qc_scorecard.py` or `run_pipeline.sh` that prints a loud,
  impossible-to-miss summary line for any `pfam_model` family flagged
  `pfam_ga_review_needed`, rather than relying on someone reading the
  full scorecard file.

## Verification (applies across phases)

`yaml.safe_load` validation (43 families, no dupes) after every
`families.yaml` edit; syntax-check touched Python/bash
(`python3 -m py_compile`); subset-test via `make_family_subset.py` + real
Docker `run_pipeline.sh all <throwaway-name>` build+benchmark before
committing any family-affecting change; compare F1/precision against the
last known-good numbers already recorded, not just "did it run." When
editing `families.yaml` or any `pipeline/*.py` file while a build is
actively running against the main directory, do the edits in a separate
`git worktree` instead (`git worktree add -b <branch> .worktrees/<name>
<base-branch>`, remove it with `git worktree remove` once merged back) —
editing a file a live pipeline run has mounted/is reading caused a real
corruption incident earlier in this project's history, and a worktree is
a physically separate checkout that avoids the problem entirely while
still letting builds run in parallel.

## Branching / PR structure

One branch + one focused PR per phase (2, 4, 6 need code changes and PRs;
1 and 3 ride along with 2's PR as follow-up commits; 5 is a decision
point, not necessarily its own PR unless it results in a real change).
File a GitHub issue for Phase 2 (new bug class, done — issue #8) and
reference it from the fixing commit, matching the pattern established in
PRs #6/#7.
