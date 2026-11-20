# Test inputs (`-profile test`)

Small inputs for a whole release, built for real in about an hour (`tests/README.md`).

- **`families.yaml`:** five families of the main `families.yaml`, copied unchanged, with the HMM
  decoys that compete with them:
  - nhaA, a plain family, and one of the real-genome check's E. coli target genes, so that check
    scores a test family;
  - betL: `decoy_from_negatives` and `betTcaiT_decoy`;
  - gor: `fadDisulfide_decoy`;
  - trkA and ktrA: the `trk_ktr_regulatory` paralog group.
- **`markers.yaml`:** three of the marker genes (ychF, pheS, ffh).
- **`pilot_strains.tsv`:** 14 of the 43 BacDive pilot strains, from 0 to 16% NaCl, most with
  complete genomes.
- **`catalog/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz`:** 6,531 genes (6.2 Mbp) of the
  GMGC v1.0 human-gut catalog (Coelho et al. 2022, Zenodo record 10473258):
  - every gene of the catalog-index slice (`docs/validation_reports/catalog_index.md`) whose best HMM match (v9) is one of the five
    families, annotated or below its cutoff (3,531, 39 of them annotated nhaA);
  - 3,000 random genes without any hit, as background.
  
  The file has the catalog's own name, because the pipeline checks it; `CATALOG_FETCH` finds it in
  place and downloads nothing.

The profile also sets:
- a random background of 500 proteins instead of 5,000;
- one BacDive chunk;

Regenerate the YAML files from the main files when those entries change. Then rerun the real test and
review the baseline (`tests/README.md`).
