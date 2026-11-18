# Catalog index: all genes or only the HMM hits?

A one-off study, not part of the pipeline. It decides whether a release's gene
catalog needs the full, competitive index of every catalog gene (`build-catalog
--index all`, about 100+ GB to build for the whole human-gut catalog and a separate
download for users), or whether the small index of the genes with an HMM hit that a
release ships (`--index hits`) is enough.

## Why the existing benchmark cannot decide it

The catalog benchmark (`docs/catalog_benchmark.md`) gave the same median F1 for both
indexes (0.892), with 170 against 270 background false positives per million reads.
But its catalog was the release's own back-translated sequences, 83% of them with an
HMM hit. A real gut catalog is about 95% genes without a hit, and those are where the
two indexes differ: with `all`, a read from such a gene maps to its own gene and is
dropped; with `hits`, it does not map and goes on to DIAMOND and the cascade. This
study tests that case directly.

## Design

1. **Genes.** A seeded random slice of the GMGC v1.0 human-gut catalog
   (`slice_genes.py`, `N_GENES` = 1 M of ~11.9 M, about 1 Gbp), in the catalog's own
   proportions.
2. **Truth.** Each gene's full-length HMM decision from `build-catalog` with the
   release's HMMs and gathering cutoffs (`catalog.annotation.tsv`): a gene annotated
   as a family is that family (on a fused gene, any family passing on it); a homolog
   below the cutoff, a decoy or a gene without a hit is no panel family.
3. **Reads** (`simulate_reads.py`). Paired 150 bp reads from the genes' nucleotides,
   started uniformly over the slice's bases, insert N(300, 30), random strand,
   substitution errors at 0.5% per base. At two divergences of the sample's strain from
   its catalog gene: 0 (the strain is in the catalog) and 2% (a close relative, below
   profile's 97% identity bar, where both indexes should fall back alike).
4. **Runs.** `osmotool profile` with its default settings on the same reads, four ways:
   `--index all`, `--index hits` (both built from the slice), the middle option below,
   and `--no_catalog`.
5. **Scores** (`score.py`). Each read's final call is rebuilt: the catalog's decision
   when it settled the read, else DIAMOND's best hit after the cascade (checked against
   the run's `gene_counts.tsv`, `verification.tsv`). Reported: false panel calls per
   million reads, split into catalog and translated calls, and recall on panel reads,
   overall (`summary.tsv`) and per family (`by_family.tsv`). The time and memory of
   every build and run are in `resources.tsv`: the full index's memory at 1 Gbp, which
   scales to the whole catalog.

Both indexes are built from the same slice, so the comparison is fair; an index of the
whole catalog is the very cost in question. A catalog call on a read from its own gene
is correct by construction (the same annotation decides both); what is measured is
where reads go and what translated search makes of those that fall through.

## Running

```bash
bash 00_setup.sh     # login node: the environment (osmotool at OSMOTOOL_REF, bwa-mem3) and release v9
bash run.sh          # submits the steps to SLURM; results in ${STUDY_DIR}/results/
bash run_neighbours.sh   # then the middle option, and the scoring of all four again
```

Settings are in `config.sh` (catalog path, slice size, reads, divergences, seed). The
release is v9 from Zenodo (`osmotool download-db`), the newest one published: the
study compares two indexes built with the same HMMs, so the release's version does not
decide it. `RELEASE=` points it at another.

## Results (2026-11-18)

Run: release v9, osmotool `a0535b4` (branch `catalog-default`), seed 20261118. The
tables are in `results/`. The per-read rebuild matched osmotool's `gene_counts.tsv` in
every cell (`verification.tsv`).

The slice: 1,001,768 genes, 755 Mbp, from a catalog of 13,788,251 genes (about
10.4 Gbp). Of the slice's genes 2,266 (0.23%) are annotated as a panel family; 22,337
are homologs below a cutoff and 8,892 decoys, so 3.3% have any HMM hit and 96.7% are
not in the hits index. 4 M reads (2 M pairs) per divergence, 12,596 of them from panel
genes.

| divergence | index | false panel calls per M reads | of which catalog / translated | recall | precision |
|---|---|---|---|---|---|
| 0 | none (`--no_catalog`) | 246.2 | 0 / 985 | 0.253 | 0.764 |
| 0 | hits | 57.2 | 65 / 164 | 0.985 | 0.982 |
| 0 | all | **0.0** | 0 / 0 | 0.983 | 1.000 |
| 2% | none | 151.2 | 0 / 605 | 0.194 | 0.802 |
| 2% | hits | 74.8 | 51 / 248 | 0.691 | 0.967 |
| 2% | all | 44.2 | 0 / 177 | 0.691 | 0.980 |

- **Either index does most of the work.** Against translated search alone, the hits
  index cuts false calls by 4x and lifts recall from 0.25 to 0.98 when the strain is
  catalogued (profile's default DIAMOND filters, 80% identity and query cover, miss most
  gut genes). Recall is the same with both indexes.
- **The full index removes the rest when the strain is catalogued**: 57 false calls per
  million reads to none. With the hits index they come from two places:
  - 65 catalog calls: reads from genes not in the hits index that align at >= 97% to a
    hit gene (cspA, trkA, opuBB, murB, mscL), the competition the full index exists for;
  - 164 translated calls, 152 of them gor, from reads of non-panel genes that fall
    through to DIAMOND. Which genes they come from is not broken down here; most likely
    gor's look-alikes (other FAD-disulfide oxidoreductases), for which v9 has no decoy.
    `fadDisulfide_decoy` (from v10 on) is built for those; if it catches them, the hits
    index's rate would drop to about 20 per million (77 of the 229). Not measured.
- **At 2% divergence** the gap narrows (75 against 44 per million). More reads settle
  than expected: 2% divergence plus 0.5% errors leaves many reads just above the 97% bar
  (2.56 M of 3.99 M non-panel reads settle with the full index).
- **Cost** (`resources.tsv`): the full index took 7.5 GB to build for 755 Mbp (about
  10 bytes per base: about 100 GB for the whole catalog) and is 2.7 GB on disk (about
  37 GB whole); profiling with it took 10.4 GB against 7.1 GB with the hits index
  (0.14 GB on disk). Each build took under 3 minutes and each profile under 2 at 32
  threads.

Not measured: a v10 release (none is published), real reads (error profiles beyond
substitutions), and the whole catalog rather than a slice.

## The middle option: the hits index plus its neighbours (2026-11-18)

The catalog false calls of the hits index come from genes it leaves out that a read can
still align to a hit gene at >= 97%. So the middle option indexes the hit genes plus the
genes without a hit that are close to one (`find_neighbours.py`): every gene without a
hit, cut into 150 bp windows every 50 bp, mapped to the hit genes with bwa-mem3 (all
placements); a gene with a window aligned at >= 90% identity over >= 45 bp is a
neighbour (looser than profile's 97% over 45 bp, for strain divergence and sequencing
errors). The neighbours have no hit, so a read they settle is dropped.

25,170 neighbours (21,344 at >= 95%, 19,592 at >= 97%; `results/neighbours.tsv.gz`):
with the 33,495 hit genes, 58,665 genes and 47 Mbp, 5.9% of the slice. Finding them
took 17 s and 9.5 GB at 32 threads; the catalog is 0.19 GB on disk (hits 0.14 GB, all
2.7 GB).

| divergence | index | false per M reads | catalog / translated | recall | precision |
|---|---|---|---|---|---|
| 0 | hits | 57.2 | 65 / 164 | 0.985 | 0.982 |
| 0 | **hits + neighbours** | **39.2** | **0** / 157 | 0.983 | 0.988 |
| 0 | all | 0.0 | 0 / 0 | 0.983 | 1.000 |
| 2% | hits | 74.8 | 51 / 248 | 0.691 | 0.967 |
| 2% | **hits + neighbours** | **61.5** | **0** / 246 | 0.691 | 0.973 |
| 2% | all | 44.2 | 0 / 177 | 0.691 | 0.980 |

- **The neighbours remove every catalog false call** (65 and 51 to none), at 1.3x the
  hits index's size instead of 19x.
- **What is left is translated, and it is gor**: 149 of the 157 at 0%, 89 of 246 at 2%.
  These reads come from genes that are not close to a gor gene in nucleotides (or they
  would be neighbours), but whose proteins DIAMOND still places on gor. The full index
  removes them only because it holds the very gene each read came from; for a strain
  that is not in the catalog it cannot. That is a decoy's job: v10's
  `fadDisulfide_decoy` targets gor's look-alikes. Without gor, the middle option leaves
  8 false calls at 0% (2 per million), against none for the full index. With v10, not
  measured.
- **At 2%** the rest is spread over gor, kdpA, galE, proP, betL and nhaA, all
  translated, so equally there with any catalog index for reads it does not settle.

Scaled to the whole catalog (13.8x the slice): about 2.5 GB on disk against about 37 GB
for the full index, and the neighbour search maps about 180 M windows, minutes at 32
threads.

## Outcome

Releases ship the middle option: osmotool `build-catalog --index neighbours`
(osmotool `f77ef36`, branch `catalog-default`), the same rule as
`find_neighbours.py`. Built with it on this slice, it found the same 25,170
neighbours, and `profile` gave gene counts identical to the study's
neighbours run. The pipeline no longer builds an index of the whole catalog.
To recheck with a v10 release once one is built: `RELEASE=<release> bash run.sh`,
then `bash run_neighbours.sh`.
