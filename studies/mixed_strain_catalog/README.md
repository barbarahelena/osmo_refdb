# Mixed-strain catalog study

One simulated sample pooling reads from genomes whose genes are in the catalog, partly in it, or absent from it, scored per group for DIAMOND alone and for the catalog path under three settle rules. It is the real-sequence counterpart of the catalog benchmark, whose back-translated reads never match real catalog genes.

```
bash run_study.sh                    # all steps; outputs in work/ (not tracked)
CATALOG_INDEX=/path/to/full_index CATALOG_TAG=gmgc_full bash run_study.sh map score verify
```

## Design

- **Genomes** (`genomes.tsv`): 12 genomes in four groups, chosen by measured catalog coverage, not by taxonomy. Panel genes of 446 BacDive genomes (all 296 gut isolates and 150 random others) were cut into 150-bp windows every 50 bp and mapped to the v11 catalog; the share of windows settling at >= 97% identity over >= 45 bp put each genome in a group: `in` (>= 0.8: E. coli K-12, Christensenella minuta, Enterococcus avium, Sutterella wadsworthensis, Campylobacter hominis), `partial` (0.2-0.8: Citrobacter europaeus, Bacteroides faecalis, Serratia nematodiphila), `out_related` (< 0.02, same genus as a catalogued strain: Enterococcus pallens, Streptococcus thoraltensis) and `out_distant` (< 0.02: Pseudomonas proteolytica, Streptomyces sulphureus). The BacDive genomes are the NCBI assemblies in `genomes.tsv` (`source`), fetched by the `genomes` step from the NCBI Datasets API; they are byte-identical to the copies the results below were computed from (BacDive phenotype validation, osmo_refdb v8). On a cluster whose compute nodes have no internet, run `bash run_study.sh genomes` on a login node first.
- **Reads**: wgsim, 150 bp pairs, 2x coverage per genome, seed = row number, default error rate 0.02 (`WGSIM_ERR`), as the release's real-genome validation. 362,292 pairs. `expected.md5` holds the checksums of the reads and the truth the results below were computed from.
- **Truth** (`truth.py`): the panel family of the gene a mate overlaps by >= 20 bp, genes from `osmotool annotate` (hmm) with the release on each genome; otherwise background. The truth is therefore the release's own full-length annotation, errors included. A wgsim name gives the fragment's outer coordinates, and either mate can sit at either end; each mate's position is found by matching its sequence to both ends.
- **Calls**: DIAMOND best hit per mate from `osmotool profile` without a catalog, with the validation's options (`--min_identity 0.40 --min_query_cover 0.50`), and the catalog placements (primary + XA) from bwa-mem3 `-h 50`. A rule's call is the catalog's decision where it settles the mate, else the DIAMOND best hit, as `compare_real_recall.py` scores the catalog path. No report levels are applied, so precision is the raw read-level precision of that scoring, not profile's reported output.
- **Methods** (`score.py`), all scored on the same reads:
  - *DIAMOND only*: no catalog.
  - *Catalogue, current rule*: osmotool 0.11.1's settle rule. A read placed on genes with different labels (a panel family and "drop") is ambiguous and gets no call. `verify` checks this rule against `osmotool profile --catalog` read by read.
  - *Catalogue, ambiguous reads to DIAMOND*: the same, but an ambiguous read is passed on to DIAMOND instead of being lost.
  - *Catalogue, ambiguous to DIAMOND + fragment rule*: also counts a homolog-labelled catalogue gene as its family when it scores like a fragment of a true gene (bits per covered HMM position at least cutoff / model length), but only when all of a read's placements then agree on one family.

  Scored families: the panel minus `profile_excluded_families.txt` (otsA, murB, mrpC, trkH, ktrB, ktrD).

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

- Strains in the catalogue: the current rule loses recall against DIAMOND only (0.54 vs 0.68) while removing 63% of the false calls. Passing ambiguous reads to DIAMOND brings recall to 0.70 with the same false calls as the current rule.
- Strains partly in the catalogue: the catalogue removes a third of the false calls without losing recall, under every rule.
- Strains not in the catalogue, related or not: the catalogue settles almost no reads (604 of ~350,000) and the results equal DIAMOND only; no reads from relatives of catalogued strains are settled as the wrong family.
- The fragment rule adds 18 correct reads for 21 false calls among catalogued strains (21 for 33 in the whole sample).
- Precision is low for every method because reads are scored on their DIAMOND best hit without report levels, as in the release's real-genome validation; profile's reported output applies report levels on top.

## Shipped index vs full index, on a SLURM cluster

The release ships the neighbours index (1,460,140 of the catalog's 13,788,251 genes: the genes with an HMM hit and the genes close to them). The full index holds every gene, so reads from genes outside the neighbours index can map to their own gene instead of to a hit gene they resemble. It is too large to build on a laptop; `hpc/` runs the study on a SLURM cluster with Apptainer, once per index, on the same reads. Both indexes are scored against v11's catalog annotation, and the full index depends only on the gene file, so the run does not wait for a new pipeline build and the index is reused by later ones.

1. Edit the first block of `hpc/env.sh`: `RELEASE` (the v11 release directory, not one a running build writes into), `CATALOG_GENES` (the gene file the release's catalog was annotated from, the build's `--catalog_genes`), `FULL_INDEX` (a new directory; the index is built there once) and `SIF_DIR`.
2. On a login node: `bash hpc/prepare.sh`. It pulls osmotool 0.11.1's image and the pipeline's wgsim image, fetches the gene file from Zenodo if it is not at `CATALOG_GENES` (one 2.99 GB member of a 19.3 GB zip, with `modules/local/catalog/fetch/templates/fetch_zenodo_member.py`) and fetches the genomes from NCBI.
3. `bash hpc/submit.sh` (partition `genoa` by default; `PARTITION` and `PARTITION_FULL` change it). It submits two jobs:
   - `study.sbatch` (16 CPUs, 32 GB, 2 h): annotation, reads, truth, DIAMOND, then the shipped index, into `work/gmgc_gut/`, with `verify` against `osmotool profile --catalog`. It checks the reads and the truth against `expected.md5`; if they differ (another wgsim build), the shipped-index numbers differ from the tables above, but both indexes are still scored on the same reads.
   - `full.sbatch` (16 CPUs, 240 GB, 12 h, after the first): the full index (`run_study.sh index`, which checks the gene file's sha256 against `catalog.json`, then runs `osmotool index-catalog`), the same reads mapped to it, scored and verified (`profile --catalog_index`), into `work/gmgc_full/`. Then `compare_indexes.py` writes `work/index_comparison.md`: recall, false calls and settled mates per strain group, for both rules and both indexes.

The full index's memory and time are not measured at this size. 240 GB and 12 h follow the pipeline's former full-index step, itself scaled from bwa-mem3 on 133 Mbp of the catalog (1.37 GB peak, 33 s). The job records what the index took (`work/gmgc_full/index.time` and the `sacct` line at the end of its log); if it runs out of memory, raise `--mem` in `full.sbatch` and submit with `PARTITION_FULL=fat_genoa`.

Checked locally with Docker: the `index` step on the repository's 6,531-gene GMGC test slice (builds the index, skips it when present, refuses a gene file with another checksum), and `truth`, `score` and `verify` with Python inside the osmotool image (`PYTHON_RUN`), which give byte-identical outputs to the run above. The Apptainer and SLURM parts have not run.

## Follow-ups

- The ambiguous-to-DIAMOND rule is in osmotool 0.12.0 (barbarahelena/osmotool#24): `profile --catalog` passes ambiguous reads to translated search, and `--catalog_only` keeps them as `catalog_ambiguous`. With that build, `profile --catalog` on this sample gives the ambiguous-to-DIAMOND column above in every group. The pipeline uses 0.12.0 since barbarahelena/osmo_refdb#35, where `catalog_benchmark.py`'s `hybrid_calls` gives ambiguous reads their translated call too. This study stays on 0.11.1: its current-rule column and `verify` step are 0.11.1's settle rule.
- The fragment rule is not implemented.
- The full index against the shipped one: see the section above; its results go here once it has run.
