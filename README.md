# osmo_refdb

osmo_refdb builds the reference database that [osmotool](https://github.com/barbarahelena/osmotool) uses to detect osmoadaptation genes, in metagenomic reads (`osmotool profile`) and in genomes (`osmotool annotate`).

One Nextflow run makes a whole release. Per gene family it fetches reference sequences from UniProt, builds an HMM and a DIAMOND database, and calibrates their cutoffs. It then decides how `osmotool profile` reports each family, benchmarks the release on held-out sequences, and checks it on real genomes. A release is a versioned folder: the databases and their cutoffs, a gene catalog for read profiling, marker genes for copies per genome, the files that tell osmotool how to report each family, and the benchmark and provenance behind them.

## Quick start

The build is a Nextflow workflow (`main.nf`; each module's `container` line and `environment.yml` give its image):

```bash
nextflow run main.nf --outdir /data/osmo_refdb/v1 --release v1    # build, benchmark, catalog and screen, BacDive, validation
nextflow run main.nf --outdir /data/osmo_refdb/v1 --release v1 -resume --bacdive_full    # the full BacDive set and the NaCl model
```

Read the catalog screen (`<outdir>/catalogs/<catalog>/screen.tsv`, printed at the end of the build) before publishing a release: a flagged family that is a real confuser class means a decoy and a rebuild under a new release name.

## What the workflow does

`main.nf` builds a whole release in one Nextflow run. The inputs are `families.yaml` (the gene panel; `--families_file` for a subset) and `markers.yaml`. The output is one release folder, `--outdir`.

1. **References.** Per family, positives and hard negatives from UniProt (plus RefSeq where a family opts in), cleaned of mislabelled sequences and length outliers, clustered at 90% identity, and split into train, calib and test.
2. **Databases.** From the train split: one HMM per family and one DIAMOND database, with decoy and fused-ORF references. The cutoffs of both are calibrated on the calib split.
3. **Benchmark.** Reads simulated from the calib split decide the cascade, the report levels, and whether the reads the gene catalog settles for each unit can be reported. Reads and proteins from the test split score the release, for `osmotool profile` and `osmotool annotate`.
4. **Gene catalog.** The catalog (GMGC human gut by default) is annotated with the release's HMMs, indexed for `osmotool profile`, and screened for families that claim genes that look like something else.
5. **Marker genes.** Ten single-copy marker genes in a DIAMOND database of their own, for copies per genome.
6. **Validation on the finished release.** The BacDive phenotype analysis (and with `--bacdive_full` the NaCl model) and the real-genome checks.
7. **Provenance.** `PROVENANCE.json` and `software_versions.yml`, last.

Nothing that is tuned sees a sequence the test benchmark scores: the HMMs and the DIAMOND database are built from train, the cutoffs, cascade and report levels are set on calib, and only the benchmark uses test.

The subworkflows are in `subworkflows/local/` and the processes in `modules/local/<tool>/<subtool>/`, each with its scripts in `templates/`. Per-process settings are in `conf/modules.config` and resources in `conf/base.config`.

Everything for a release is written to `--outdir`, one folder per release, so you can keep multiple versions around and reproduce any of them independently.

### Options and outputs

Every build option, with its default, is in [docs/build_options.md](docs/build_options.md). What each process publishes into the release folder is in [docs/outputs.md](docs/outputs.md).

## Gene families

The panel, every `families.yaml` field (including `report_level` and the reason codes in the report-level manifests), and how to add and test a new family are in [docs/FAMILIES_SUMMARY.md](docs/FAMILIES_SUMMARY.md).

## Validation

Every release is checked in six ways. The first five run in every build; the taxonomy split is a separate build. Where each check's results are published is in [docs/outputs.md](docs/outputs.md).

- **Held-out benchmark (simulated reads).** Reads simulated from the test split, which no step that is tuned has seen, plus hard negatives from close paralogs and a random bacterial background. `osmotool profile` is scored per family and per reporting unit, for DIAMOND alone, DIAMOND with the HMM cascade, and HMM alone. A real metagenome can be added as extra background (`--background_r1/_r2`).
- **Annotate benchmark (full-length proteins).** `osmotool annotate` on the test split's proteins, their hard negatives and the background, per family and paralog group, for HMM and DIAMOND.
- **Catalog benchmark.** `osmotool profile` with a gene catalog on the same test reads, once with every read's own gene catalogued (the best case) and once without it (reads from strains the catalog does not hold).
- **Real genomes.** On published reference genomes (`--reference_genomes_dir`): whether osmotool calls the multi-subunit systems that are known to be there and not those known to be absent. On genomes with curated target genes: read-level recall and precision on reads simulated from the genome, with and without the catalog.
- **Phenotype (BacDive).** Whether the panel's genes predict salt tolerance across BacDive strains better than random gene panels of the same size, with a per-family ranking, cross-validation and a within-domain control for salt-in strategists. With `--bacdive_full` the same data train the NaCl model ([docs/nacl_model.md](docs/nacl_model.md)).
- **Unseen genera (`--split_mode taxonomy`).** A build that holds out whole genera instead of random sequences, to test detection of divergent homologs.

## Licence and citations

The code in this repository is MIT-licensed (`LICENSE`).

A release (the DIAMOND database, HMMs, reference sequences, gene catalog and results that `osmotool` downloads) is not code, and is licensed separately, under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Reuse needs attribution to the sources it is built from:

- **UniProt** for the reference sequences, and through them the HMMs and the DIAMOND database: The UniProt Consortium, *UniProt: the Universal Protein Knowledgebase in 2025*, Nucleic Acids Res 53(D1):D609-D617 (2025), [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010).
- **GMGC** for the gene catalog under `catalogs/`: Coelho et al., *Towards the biogeography of prokaryotic genes*, Nature 601:252-256 (2022), [doi:10.1038/s41586-021-04233-4](https://doi.org/10.1038/s41586-021-04233-4).
- **The prevalence-filtered catalogs** the gene set is taken from: Ducarmon et al., *Cayman enables large-scale analysis of gut microbiome carbohydrate-active enzyme repertoires*, Nat Microbiol 11:1739-1753 (2026), [doi:10.1038/s41564-026-02318-2](https://doi.org/10.1038/s41564-026-02318-2); data at [zenodo.org/records/10473258](https://zenodo.org/records/10473258).

Each release's `PROVENANCE.json` records both licences and these sources under `licence`, so they travel with the data.

`CITATIONS.md` lists the tools and databases the pipeline uses, with the references to cite for them.
