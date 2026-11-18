# Real-genome validation

The release benchmark measures DIAMOND vs HMM performance on reads simulated from UniProt-derived
sequences. This is a separate, complementary check: does the database call genes correctly on real,
independently annotated genomes, judged against annotation this project did not produce (NCBI RefSeq,
UniProt)?

It runs for every release, as part of the Nextflow workflow (`subworkflows/local/real_genome_validation`),
after the gate. Results go to `releases/<release>/real_genome_validation/`.

## What runs

1. **System calls on the nine reference genomes** (osmotool's `test_genomes/`, accessions in its
   README; pass the directory as `reference_genomes_dir`). `osmotool annotate` on each genome, then
   `check_reference_systems.py` against `assets/real_genome_validation/expected_systems.tsv`: present,
   absent or `info` (reported, not scored), with the basis for each. A failed expectation is written
   to `reference_systems_check.txt`; it does not stop the workflow.
2. **Read-level recall and precision** on the genomes with curated target genes
   (`assets/real_genome_validation/target_orfs.tsv`; the genomes are `<test_genomes_dir>/<genome>.fna`,
   default `assets/test_genomes/` in the repo). Reads are simulated from the genome with `wgsim`, whose headers
   carry the true coordinates, so the ground truth is independent of the database. DIAMOND
   (`osmotool profile --no_catalog`) and the HMMs (`orfm` + `hmmscan`) are scored per family against it
   (`real_read_truth.py`, `compare_real_recall.py`) in `<genome>.real_recall_comparison.tsv`. So is the
   `catalog` method: `osmotool profile` with the release's gene catalog and the index it ships, as
   profile runs by default. A read the catalog settles takes its decision, every other read its
   DIAMOND call. These are the only real-genome reads the catalog is tested on; a genome whose strain is
   not in the catalog (likely for the halophiles in a gut catalog) shows what the fallback gives.

## The manual part, done once

Deciding which ORF is which gene is manual on purpose: it is what keeps the check independent of the
tool. It is done once per genome, not per release:

1. Annotate the genome with any release (`osmotool annotate ... --keep_proteins`).
2. Cross-reference each hit against independent annotation (RefSeq gene names, UniProt EC numbers,
   review status) to confirm which calls are real; leave out known false positives.
3. Add one row per confirmed gene to `assets/real_genome_validation/target_orfs.tsv`
   (`genome`, `family`, Prodigal ORF id from the `.prodigal.faa` header) and put the genome in
   `assets/test_genomes/<genome>.fna`.

Every later release is then validated against the same curated truth. Only the source genome FASTAs
belong in `assets/test_genomes/`; the outputs are git-ignored release files.

The read-recall and system-check scripts can also be run by hand (`python real_read_truth.py --help`, and so on); the target coordinates are worked out in `VALIDATION_TARGETCOORDS` itself.
