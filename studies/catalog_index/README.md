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
4. **Runs.** `osmotool profile` with its default settings on the same reads, three ways:
   `--index all`, `--index hits` (both built from the slice) and `--no_catalog`.
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
```

Settings are in `config.sh` (catalog path, slice size, reads, divergences, seed). The
release is v9 from Zenodo (`osmotool download-db`), the newest one published: the
study compares two indexes built with the same HMMs, so the release's version does not
decide it. `RELEASE=` points it at another.

## Results

To be filled in from `${STUDY_DIR}/results/` once the run is done.
