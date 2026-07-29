# families.yaml Summary

`families.yaml` is the single source of truth for osmo_refdb's gene family panel — 43 families as of 2026-07-24. Adding a family here and rerunning `bash run_pipeline.sh build` is the only change needed; no other files require edits.

## Schema

| Field | Purpose |
|---|---|
| `name` | Short gene/family code — FASTA tag and HMM model name (unique) |
| `description` | Free text for manifests/reports |
| `positive_query` | UniProt REST query for the true-positive set |
| `negative_query` | UniProt REST query for a "hard negative" set (prefer `xref:pfam-PFxxxxx` over free-text `family:`, which silently zero-hits) |
| `negative_pfam` | Pfam accession documented for the negative query |
| `scope: annotate_only` | Excludes a family from `osmotool profile` output while still building/searching it (e.g. housekeeping genes) |
| `decoy_from_negatives` | Adds the QC'd negative set to the DIAMOND db as a `<family>_decoy`, letting a confusable named paralog win reads away from this family instead of relying on a score cutoff |
| `trim_gt` | Per-family override of trimAl's gap threshold (default 0.8) — lower values keep more variable/flanking columns |
| `pfam_model` | Use Pfam's own curated HMM directly instead of building one from our fetched sequences (only when a per-gene InterPro check shows the Pfam family is gene-specific, not a broad fold) |
| `fusion_partner` / `fusion_marker_pfam` | Declares that this gene occurs as a single fused ORF with a named partner family in some lineages, plus the marker domain that confirms a real fusion vs. a length outlier. Also drives a negative-side pre-filter (a fused ORF under a bare locus tag matches this family's own negative_query just as easily as a true hard negative) and query-time exclusion — see mrpA/mrpB, otsA/otsB below |
| `max_positive_override` | Caps the positive fetch size (used for cspA, ~80k members in Pfam) |
| `max_negative_override` | Raises the negative fetch size above the default 1000 cap (only mazG so far — see below) |
| `negative_median_override` | Pins the negative-set length-filter median to a fixed value instead of recomputing it from whatever the (capped) fetch contains this run (only mazG so far — see below) |

## Family groups

**Ectoine synthesis** — ectA, ectB, ectC (steps 1–3; all carry numbered-paralog synonyms)

**Trehalose synthesis** — otsA, otsB (confirmed fused-ORF variant, same pattern as mrpA/mrpB — see below)

**Glutathione redox cycle** — gshA, gshB, gshF (bifunctional gshA/gshB fusion-equivalent found in Streptococcus/Enterococcus/Listeria/Clostridium), gor

**Choline → glycine betaine oxidation** — betA, betB (Ng et al. 2023 mBio, in vivo gut validation)

**Compatible-solute transporters**
- betL (BCCT family; folds betL/betP/opuD/betS orthologs into one positive class, unlike the paralog groups below)
- proX, opuAC, opuBC, opuCC (ABC substrate-binding subunits — real paralogs, kept as 4 separate families that cross-exclude each other's negatives)
- proP (MFS transporter)
- opuAA, opuBA, opuCA (ABC ATPase subunits); opuAB, opuBB, opuCB (ABC permease subunits)

**Mrp/Mnh Na+/H+ antiporter complex (Bacillota nhaA substitute)** — mrpA/mnhA, mrpB/mnhB, mrpC/mnhC, mrpD/mnhD, mrpE/mnhE, mrpF/mnhF/PhaF, mrpG/mnhG/yufB (mrpA+mrpB have a confirmed fused-ORF variant)

**Other Na+/H+ and K+ transport** — nhaA (Proteobacteria-type antiporter); kdpA (inducible high-affinity K+ uptake); trkA/trkH (Proteobacteria-named constitutive K+ uptake) and ktrA/ktrB/ktrD (Firmicutes-named counterpart — kept separate due to a hard HMM tie-breaking bug when Pfam-sharing siblings adopt the same model)

**Mechanosensitive channels** — mscL, mscS

**Housekeeping / co-occurrence markers (Culligan et al. 2012 salt-tolerance locus)** — galE, mazG (only family in the panel needing `max_negative_override`/`negative_median_override` — see below), murB (`scope: annotate_only`)

**Cold shock** — cspA (the only family anchored on a bare Pfam accession rather than a gene symbol, to capture all csp paralogs per genome)

## Recurring design patterns worth knowing

- **Numbered-paralog gap (issue #8):** most families also match `<gene>1`/`<gene>2` symbols — real, reviewed UniProt gene names for organisms carrying two distinct paralogs, not noise.
- **Ortholog vs. paralog handling:** cross-organism orthologs under different names (betL/betP/opuD/betS, mrpA/mnhA) are folded into one positive class; same-organism paralogs with distinct specificities (proX vs. opuAC/opuBC/opuCC, ktrB vs. ktrD) are kept as separate families that exclude each other from their own negative pools.
- **Promiscuous-domain caution:** several families (ectB, mscS, galE, mrpC, betB, gshA) are flagged to check `qc_scorecard.tsv`/`cutoff_manifest.tsv` before trusting the HMM, since their negative-anchoring Pfam domain is broader than the gene itself.
- **`pfam_model` is applied selectively**, not panel-wide — only after confirming via InterPro (low domain-architecture count, name matches the gene 1:1) that Pfam's own model is gene-specific, and never for a Pfam family shared by more than one distinct sibling family in this panel (breaks HMM discrimination — see trkH/ktrB/ktrD).
- **`decoy_from_negatives`** is currently only set for betL (betT/caiT), where calibration showed a single score cutoff structurally cannot separate the two.
- **A fusion partner can contaminate a family's *negative* pool, not just its positive fetch.** A genuinely fused ORF under a bare locus tag (no recognizable gene symbol) matches a family's own `negative_query` just as easily as a true hard negative, since it carries the domain the query is anchored on. Confirmed on both declared fusion pairs: mrpB (~68% of its negative pool) and otsA/otsB (~41–59%) — in both cases the contamination was the *majority* of the pool, which silently became the "normal" length center and caused the median-based length-outlier filter to discard the true, correctly-sized hard negatives instead of the contamination. Checked mrpF/mrpG/gshA/gshB/trkH/ktrB/ktrD for the same pattern (2026-07-29) — none showed it; their contamination is the plain "unlabeled real ortholog" kind below, not a hidden fusion.
- **Negative-set length-filter median instability (capped-fetch sample not representative of the true population) is real but rare.** Checked systematically across every family whose negative fetch hits the default cap (38 of 43) — only mazG showed it (215aa at the default n=1000 cap vs. 116aa at n=3000); every other family's median moved ~2% or less regardless of fetch size. Not treated as a general mechanism — just a documented one-off fix on mazG's own entry.

## Deferred families (not added, 2026-07-16)

proV (PF00005), proW (PF00528), gbsA (PF00171), gbsB (PF00465) — all anchor on Pfam domains too broad (shared by huge numbers of unrelated bacterial proteins) for reliable single-domain HMM/DIAMOND detection. Revisiting these needs hand-picked negatives rather than a raw Pfam pool.
