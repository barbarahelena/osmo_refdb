# Real-genome validation

`run_pipeline.sh`'s `benchmark` step measures DIAMOND vs HMM performance on
reads simulated from UniProt-derived sequences (synthetic, though built
from real protein sequences). This directory holds a separate, complementary
check: does the reference database actually call genes correctly on a real,
independently-annotated genome, cross-referenced against annotation this
project didn't produce (NCBI RefSeq / UniProt), not just against itself?

Unlike the main pipeline, this workflow is **semi-manual by design**: step 2
below requires a human to look at real annotation and decide which ORF is
which gene. That's the point — it's an independent check, not something
that can be automated without losing the independence.

## Workflow

1. **Get a real, annotated genome.** Pick something with an existing NCBI
   RefSeq or UniProt annotation to cross-reference against. Past runs used
   E. coli K-12 (`GCF_000005845.2`/`NC_000913.3`), *Rhodopirellula baltica*
   (`NC_005027`), and *Haloferax volcanii* (`GCF_000025685.1`, an archaeon —
   zero training-set overlap with the mostly-bacterial reference sequences).
   Save the genome FASTA to `test_genomes/<organism>.fna`.

2. **Annotate it and manually confirm the calls.**
   ```bash
   osmotool annotate releases/<release>/osmo_refdb.dmnd test_genomes/<organism>.fna \
       --method both --hmm_db releases/<release>/hmms/osmo_refdb.hmm \
       --keep_proteins --keep_aln --out_prefix releases/<release>/real_genome_validation/<organism>/<organism>
   ```
   Look at the resulting `.hmmscan.tblout` / `.blastp.tsv` calls and
   cross-reference each hit's protein against independent annotation (NCBI
   RefSeq gene names, UniProt EC numbers, review status) to confirm which
   calls are real. This is the actual verification step — the tool's own
   output is never trusted circularly as its own ground truth.

3. **Build a coordinates table for the confirmed genes.** Edit
   `extract_target_coords.py`'s `TARGET_ORFS` dict to map each family to
   the Prodigal ORF ID you confirmed in step 2 (from the retained
   `--keep_proteins` protein FASTA's headers), then run:
   ```bash
   python extract_target_coords.py \
       --prodigal-faa releases/<release>/real_genome_validation/<organism>/<organism>.prodigal.faa \
       --out releases/<release>/real_genome_validation/<organism>/gene_regions.tsv
   ```

4. **Simulate short reads directly from the genome** with `wgsim` (its
   embedded read-coordinate headers give exact ground truth with no
   k-mer matching needed, unlike the main pipeline's
   `09b_compute_read_truth.py`, which has to support simulators that don't
   expose coordinates):
   ```bash
   wgsim -N 20000 -1 150 -2 150 test_genomes/<organism>.fna \
       releases/<release>/real_genome_validation/<organism>/reads_R1.fastq \
       releases/<release>/real_genome_validation/<organism>/reads_R2.fastq
   ```

5. **Run DIAMOND (`osmotool profile`) and HMM (`orfm` + `hmmscan`) on those
   reads**, then score both against coordinate-based ground truth:
   ```bash
   python real_read_truth.py \
       --r1 releases/<release>/real_genome_validation/<organism>/reads_R1.fastq \
       --r2 releases/<release>/real_genome_validation/<organism>/reads_R2.fastq \
       --regions releases/<release>/real_genome_validation/<organism>/gene_regions.tsv \
       --read-length 150 \
       --out releases/<release>/real_genome_validation/<organism>/read_truth.tsv

   python compare_real_recall.py \
       --truth releases/<release>/real_genome_validation/<organism>/read_truth.tsv \
       --diamond-blastx releases/<release>/real_genome_validation/<organism>/reads.blastx.tsv \
       --hmm-tblout releases/<release>/real_genome_validation/<organism>/reads.hmmscan.tblout \
       --out releases/<release>/real_genome_validation/<organism>/real_recall_comparison.tsv
   ```

Write everything under `releases/<release>/real_genome_validation/<organism>/`
(git-ignored, like other release output) rather than committing one-off
intermediate files here — only the source genome FASTAs belong in
`test_genomes/`.
