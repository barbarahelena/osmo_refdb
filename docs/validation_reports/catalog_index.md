# Catalog index: which genes the release's catalog indexes

Tested on 2026-11-18. It decided the index a release's gene catalog ships with:
`osmotool build-catalog --index neighbours`.

## The question

`osmotool profile` with a catalog settles a read when it aligns at ≥97% identity over ≥45 bp to a
catalog gene; every other read goes to DIAMOND and the cascade. Which genes the index holds decides
where a read from a gene *without* a panel HMM hit goes:

- **All genes** (`--index all`): the read maps to its own gene and is dropped. But the index of the
  whole human-gut catalog takes about 100 GB to build and would be a separate download for users.
- **Only the genes with an HMM hit** (`--index hits`), sub-cutoff homologs and decoys included: small
  enough to ship. But a read from any other gene either falls through to DIAMOND or settles on a
  hit gene it is close to.

The catalog benchmark (`catalog_benchmark.md`) could not decide it. Its catalog was the release's own
back-translated sequences, 83% of them with a hit. In a real gut catalog about 97% of the genes have
no hit, and those are exactly the genes where the two indexes differ.

## What was tested

**Genes:**
- A seeded random slice of the GMGC v1.0 human-gut catalog: 1,001,768 of its 13,788,251 genes
  (755 Mbp of about 10.4 Gbp), in the catalog's own proportions.
- Of these, 2,266 (0.23%) are annotated as a panel family. Another 22,337 are homologs below a
  cutoff and 8,892 are decoys, so 3.3% have any hit and 96.7% are not in the hits index.

**Truth:**
- Each gene's full-length HMM decision from `build-catalog` with release v9's HMMs and gathering
  cutoffs. A gene annotated as a family is that family; on a fused gene, any family passing on it.
- Every other gene (a homolog below the cutoff, a decoy, or no hit at all) is no panel family.
- v9 was the newest published release. All indexes were built with the same HMMs, so the release's
  version doesn't decide the comparison.

**Reads:**
- Paired 150 bp reads simulated from the genes' nucleotides, not from back-translated proteins:
  started uniformly over the slice's bases, insert N(300, 30), random strand, 0.5% substitution
  errors.
- 2 million pairs at each of two divergences of the sample's strain from its catalog gene:
  - 0%: the strain is in the catalog;
  - 2%: a close relative.
- 12,596 of the 4 million reads at each divergence come from panel genes.

**Four ways of profiling the same reads,** with `osmotool profile`'s default settings:
1. no catalog (`--no_catalog`);
2. the hits index;
3. the hits index plus the hit genes' **neighbours**: every gene without a hit that is close to a
   hit gene, found by cutting it into 150 bp windows every 50 bp, mapping them to the hit genes with
   bwa-mem3 (all placements), and keeping it if a window aligns at ≥90% over ≥45 bp. That's looser
   than profile's 97%, to allow for strain divergence and sequencing errors. A read that settles on
   a neighbour is dropped, since it has no hit;
4. the index of all genes.

**Scoring:**
- Each read's final call was rebuilt: the catalog's decision when the catalog settled the read,
  otherwise DIAMOND's best hit after the cascade. It was checked against osmotool's own
  `gene_counts.tsv`, and matched in every cell.
- A false panel call is a call to a family the read's gene is not.
- Recall is over the reads from panel genes.

osmotool `a0535b4` (and `f77ef36` for the neighbours index), seed 20261118.

## Results

| divergence | index | false panel calls per M reads | of which catalog / DIAMOND calls | recall | precision |
|---|---|---|---|---|---|
| 0% | no catalog | 246.2 | 0 / 985 | 0.253 | 0.764 |
| 0% | hits | 57.2 | 65 / 164 | 0.985 | 0.982 |
| 0% | **hits + neighbours** | **39.2** | **0** / 157 | 0.983 | 0.988 |
| 0% | all | 0.0 | 0 / 0 | 0.983 | 1.000 |
| 2% | no catalog | 151.2 | 0 / 605 | 0.194 | 0.802 |
| 2% | hits | 74.8 | 51 / 248 | 0.691 | 0.967 |
| 2% | **hits + neighbours** | **61.5** | **0** / 246 | 0.691 | 0.973 |
| 2% | all | 44.2 | 0 / 177 | 0.691 | 0.980 |

(Counts are reads out of 4 million.)

- **Any catalog index does most of the work.**
  - Against translated search alone, it cuts false calls by four and lifts recall from 0.25 to 0.98
    when the strain is catalogued. Profile's default DIAMOND filters, 80% identity and query cover,
    miss most gut genes.
  - Recall is the same with every index.
- **The hits index's own catalog false calls come from missing competitors.** 65 at 0% and 51 at 2%
  are reads from genes left out of the index that still align at ≥97% to a hit gene, spread over
  cspA, trkA, opuBB, murB and mscL.
- **Adding the 25,170 neighbours removes every one of them.** The index then holds 58,665 genes,
  1.3× the hits index, against 19× for the index of all genes.
- **What remains is DIAMOND calls, and mostly gor:** 149 of the 157 at 0% and 89 of the 246 at 2%.
  - These reads come from genes that aren't close to a gor gene in nucleotides, or they would be
    neighbours, but whose proteins DIAMOND still places on gor.
  - The index of all genes removes them only because it holds the very gene each read came from. For
    a strain that isn't in the catalog it can't, so a decoy has to handle them either way.
  - Where they come from wasn't broken down. Most likely they're gor's look-alikes (other
    FAD-disulfide oxidoreductases), which v9 has no decoy for; v10's `fadDisulfide_decoy` targets
    them.
  - Without gor, the neighbours index leaves 8 false calls at 0% (2 per million), against none for
    the index of all genes.
- **At 2% divergence**, more reads settle than the 97% bar suggests: 2% divergence plus 0.5% errors
  leaves many reads just above it. The DIAMOND false calls left are spread over gor, kdpA, galE, proP,
  betL and nhaA.

**Size and cost on the slice,** scaled linearly to the whole catalog (13.8× the slice):

| index | genes | on disk | whole catalog, scaled | peak memory to build |
|---|---|---|---|---|
| hits | 33,495 | 0.14 GB | about 1.9 GB | 0.7 GB |
| hits + neighbours | 58,665 | 0.19 GB | about 2.5 GB | 9.5 GB (the window mapping) |
| all | 1,001,768 | 2.7 GB | about 37 GB | 7.5 GB (about 10 bytes per base: about 100 GB whole) |

Each build took under 3 minutes and each profile under 2, at 32 threads. Profiling took 10.4 GB with
the index of all genes against 7.1–7.2 GB with the others.

## Outcome

- **Releases ship the neighbours index:** osmotool `build-catalog --index neighbours`, the same rule
  as tested. Rebuilt with it on the same slice, it found the same 25,170 neighbours, and `profile`
  gave identical gene counts.
- **No whole-catalog index:** the pipeline no longer builds one.

**Not measured:**
- a v10 release, which would show whether `fadDisulfide_decoy` removes the gor calls;
- real reads, with error profiles beyond substitutions;
- the whole catalog rather than a slice.
