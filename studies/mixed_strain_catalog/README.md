# Mixed-strain catalog study

One simulated sample pooling reads from genomes whose genes are in the catalog, partly
in it, or absent from it, scored per group for DIAMOND alone and for the catalog path
under three settle rules. It is the real-sequence counterpart of the catalog benchmark,
whose back-translated reads never match real catalog genes.

```
bash run_study.sh                    # all steps; outputs in work/ (not tracked)
CATALOG_INDEX=/path/to/full_index CATALOG_TAG=gmgc_full bash run_study.sh map score verify
```

## Design

- **Genomes** (`genomes.tsv`): 12 genomes in four groups, chosen by measured catalog
  coverage, not by taxonomy. Panel genes of 446 BacDive genomes (all 296 gut isolates and
  150 random others) were cut into 150-bp windows every 50 bp and mapped to the v11
  catalog; the share of windows settling at >= 97% identity over >= 45 bp put each
  genome in a group: `in` (>= 0.8: E. coli K-12, Christensenella minuta, Enterococcus
  avium, Sutterella wadsworthensis, Campylobacter hominis), `partial` (0.2-0.8:
  Citrobacter europaeus, Bacteroides faecalis, Serratia nematodiphila), `out_related`
  (< 0.02, same genus as a catalogued strain: Enterococcus pallens, Streptococcus
  thoraltensis) and `out_distant` (< 0.02: Pseudomonas proteolytica, Streptomyces
  sulphureus).
  The BacDive genomes are the NCBI assemblies in `genomes.tsv` (`source`), fetched by the
  `genomes` step from the NCBI Datasets API; they are byte-identical to the copies the
  results below were computed from (BacDive phenotype validation, osmo_refdb v8). On a
  cluster whose compute nodes have no internet, run `bash run_study.sh genomes` on a
  login node first.
- **Reads**: wgsim, 150 bp pairs, 2x coverage per genome, seed = row number, default
  error rate 0.02 (`WGSIM_ERR`), as the release's real-genome validation. 362,292 pairs.
- **Truth** (`truth.py`): the panel family of the gene a mate overlaps by >= 20 bp, genes
  from `osmotool annotate` (hmm) with the release on each genome; otherwise background.
  The truth is therefore the release's own full-length annotation, errors included. A
  wgsim name gives the fragment's outer coordinates, and either mate can sit at either
  end; each mate's position is found by matching its sequence to both ends.
- **Calls**: DIAMOND best hit per mate from `osmotool profile` without a catalog, with
  the validation's options (`--min_identity 0.40 --min_query_cover 0.50`), and the
  catalog placements (primary + XA) from bwa-mem3 `-h 50`. A rule's call is the
  catalog's decision where it settles the mate, else the DIAMOND best hit, as
  `compare_real_recall.py` scores the catalog path. No report levels are applied, so
  precision is the raw read-level precision of that scoring, not profile's reported
  output.
- **Methods** (`score.py`), all scored on the same reads:
  - *DIAMOND only*: no catalog.
  - *Catalogue, current rule*: osmotool 0.11.1's settle rule. A read placed on genes
    with different labels (a panel family and "drop") is ambiguous and gets no call.
    `verify` checks this rule against `osmotool profile --catalog` read by read.
  - *Catalogue, ambiguous reads to DIAMOND*: the same, but an ambiguous read is passed on
    to DIAMOND instead of being lost.
  - *Catalogue, ambiguous to DIAMOND + fragment rule*: also counts a homolog-labelled
    catalogue gene as its family when it scores like a fragment of a true gene (bits per
    covered HMM position at least cutoff / model length), but only when all of a read's
    placements then agree on one family.

  Scored families: the panel minus `profile_excluded_families.txt` (otsA, murB, mrpC,
  trkH, ktrB, ktrD).

## Results (v11, shipped neighbours index, 2026-11-26)

**Reads from panel genes that get the right family (recall)**

| Strains | DIAMOND only (no catalogue) | Catalogue, current rule | Catalogue, ambiguous reads to DIAMOND | Catalogue, ambiguous to DIAMOND + fragment rule |
|---|---:|---:|---:|---:|
| Strains in the catalogue (5 genomes) | 0.68 (709/1051) | 0.54 (567/1051) | 0.70 (734/1051) | 0.72 (752/1051) |
| Strains partly in the catalogue (3) | 0.73 (756/1030) | 0.74 (766/1030) | 0.76 (779/1030) | 0.76 (782/1030) |
| Not in the catalogue, same genus as a catalogued strain (2) | 0.61 (318/521) | 0.61 (318/521) | 0.61 (318/521) | 0.61 (318/521) |
| Not in the catalogue, distant (2) | 0.76 (885/1165) | 0.76 (885/1165) | 0.76 (885/1165) | 0.76 (885/1165) |
| Whole sample (12) | 0.71 (2668/3767) | 0.67 (2536/3767) | 0.72 (2716/3767) | 0.73 (2737/3767) |

**False calls (reads given a family they are not from)**

| Strains | DIAMOND only (no catalogue) | Catalogue, current rule | Catalogue, ambiguous reads to DIAMOND | Catalogue, ambiguous to DIAMOND + fragment rule |
|---|---:|---:|---:|---:|
| Strains in the catalogue (5 genomes) | 1,685 | 619 | 623 | 644 |
| Strains partly in the catalogue (3) | 1,963 | 1,340 | 1,340 | 1,352 |
| Not in the catalogue, same genus as a catalogued strain (2) | 769 | 768 | 768 | 768 |
| Not in the catalogue, distant (2) | 1,864 | 1,856 | 1,856 | 1,856 |
| Whole sample (12) | 6,281 | 4,583 | 4,587 | 4,620 |

- Strains in the catalogue: the current rule loses recall against DIAMOND only (0.54
  vs 0.68) while removing 63% of the false calls. Passing ambiguous reads to DIAMOND
  brings recall to 0.70 with the same false calls as the current rule.
- Strains partly in the catalogue: the catalogue removes a third of the false calls
  without losing recall, under every rule.
- Strains not in the catalogue, related or not: the catalogue settles almost no reads
  (604 of ~350,000) and the results equal DIAMOND only; no reads from relatives of
  catalogued strains are settled as the wrong family.
- The fragment rule adds 18 correct reads for 21 false calls among catalogued strains
  (21 for 33 in the whole sample).
- Precision is low for every method because reads are scored on their DIAMOND best hit
  without report levels, as in the release's real-genome validation; profile's reported
  output applies report levels on top.

## Follow-ups

- The ambiguous-to-DIAMOND rule is implemented in osmotool (barbarahelena/osmotool#24,
  for 0.12.0): `profile --catalog` passes ambiguous reads to translated search, and
  `--catalog_only` keeps them as `catalog_ambiguous`. With that build, `profile --catalog`
  on this sample gives the ambiguous-to-DIAMOND column above in every group. When the
  build moves to that image, `catalog_benchmark.py`'s `hybrid_calls` must give ambiguous
  reads their translated call too.
- The fragment rule is not implemented.
- The full GMGC catalogue (13.8 M genes) against the shipped neighbours index: on the
  HPC, with `CATALOG_INDEX` set to the full index (second command at the top).
