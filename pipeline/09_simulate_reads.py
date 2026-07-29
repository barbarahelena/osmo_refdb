#!/usr/bin/env python3
"""
09_simulate_reads.py — simulate short reads from held-out test positives +
hard negatives, for benchmarking DIAMOND vs HMM.

Steps:
  1. Load each family's held-out positive test set (refs/<family>.positive.test.faa,
     produced upstream by 03_split_train_test.py so it was never used to
     build the HMM or DIAMOND db) and hard-negative set (refs/<family>.negative.faa).
  2. Optionally mix in real background metagenome reads (--background) that
     contain no known osmoadaptation genes, to estimate a realistic
     false-positive rate in community context.
  3. Simulate paired-end Illumina-like reads from the held-out protein
     sequences (embedded in random flanking nucleotide context so ORFs sit
     mid-read, mimicking real gene fragments in metagenomic reads) using
     InSilicoSeq if available, else wgsim.
  4. Write reads + a truth table (read_id -> family/negative/background,
     source UniProt ID) for later scoring by 11_compute_metrics.py.

Reproducibility: wgsim (-S) and InSilicoSeq (--seed) both need an explicit
seed passed at the CLI to be deterministic -- neither seeds itself by
default, so two runs against identical held-out input sequences previously
produced slightly different actual reads (different read counts, sequencing
errors, sampled positions) purely from the simulator's own unseeded
internal RNG. Each family/label's simulator invocation now gets a
deterministic seed drawn from this script's own RANDOM_SEED-seeded rng, so
a re-run against unchanged upstream data reproduces byte-identical reads.

Usage:
  python 09_simulate_reads.py --refs refs --families families.yaml \\
      --reads-per-sequence 20 --out results/reads
"""

from __future__ import annotations

import argparse
import csv
import random
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

RANDOM_SEED = 42
FLANK_LEN = 300  # bp of random flanking DNA on each side of the CDS


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def load_decoy_family_names(path: Path) -> set[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return {fam["name"] for fam in data["families"] if fam.get("decoy_from_negatives")}


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header, sequence) tuples."""
    records = []
    header, seq = None, []
    if not path.exists():
        return records
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def protein_to_dna(protein: str) -> str:
    """Reverse-translate a protein to a DNA sequence using a simple codon table
    (most-frequent codon per amino acid; sufficient for simulated reads, not
    intended to be biologically representative of real codon usage)."""
    codon_table = {
        "A": "GCC", "R": "CGC", "N": "AAC", "D": "GAC", "C": "TGC",
        "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
        "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCG",
        "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
        "X": "NNN", "*": "",
    }
    return "".join(codon_table.get(aa, "NNN") for aa in protein.upper())


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_genomic_context(protein: str, rng: random.Random) -> str:
    cds = protein_to_dna(protein)
    return random_dna(FLANK_LEN, rng) + cds + random_dna(FLANK_LEN, rng)


def which_simulator() -> str:
    if shutil.which("iss"):
        return "iss"
    if shutil.which("wgsim"):
        return "wgsim"
    print("ERROR: neither InSilicoSeq ('iss') nor 'wgsim' found on PATH.",
          file=sys.stderr)
    sys.exit(1)


def simulate_with_iss(genome_fasta: Path, out_prefix: Path, n_reads: int, seed: int) -> None:
    subprocess.run(
        ["iss", "generate", "--genomes", str(genome_fasta),
         "--n_reads", str(n_reads), "--model", "hiseq",
         "--seed", str(seed),
         "--output", str(out_prefix)],
        check=True,
    )


def simulate_with_wgsim(genome_fasta: Path, out_prefix: Path, n_reads: int,
                         seed: int, read_length: int = 150) -> None:
    # wgsim needs per-sequence coverage; simplest robust approach is to
    # concatenate all contigs and let wgsim sample reads genome-wide.
    r1 = f"{out_prefix}_R1.fastq"
    r2 = f"{out_prefix}_R2.fastq"
    n_pairs = max(1, n_reads // 2)
    subprocess.run(
        ["wgsim", "-N", str(n_pairs), "-1", str(read_length), "-2", str(read_length),
         "-S", str(seed),
         str(genome_fasta), r1, r2],
        check=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs", required=True, type=Path, help="Path to osmo_refdb/refs")
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--reads-per-sequence", type=int, default=20,
                     help="Simulated read pairs per held-out sequence")
    ap.add_argument("--read-lengths", type=str, default=None,
                     help="Comma-separated read lengths in bp (e.g. "
                          "'100,150,250,300') to simulate at, for "
                          "read-length-stratified metrics. Forces wgsim "
                          "(exact length control) regardless of whether "
                          "InSilicoSeq is available. Default: single run at "
                          "whichever simulator's normal default length.")
    ap.add_argument("--background", nargs=2, metavar=("R1", "R2"), default=None,
                     help="Optional real paired-end background FASTQ files "
                          "with no known osmoadaptation content -- run "
                          "through the same DIAMOND/HMM benchmark as every "
                          "other sample, scored as pure false-positive risk "
                          "in background_fpr.tsv (11_compute_metrics.py), "
                          "since they don't belong to any gene family")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rng = random.Random(RANDOM_SEED)
    families = load_family_names(args.families)
    decoy_families = load_decoy_family_names(args.families)
    args.out.mkdir(parents=True, exist_ok=True)
    read_lengths = ([int(x) for x in args.read_lengths.split(",")]
                     if args.read_lengths else None)
    genome_dir = args.out / "genomic_context"
    genome_dir.mkdir(exist_ok=True)

    simulator = which_simulator()

    truth_rows = []

    for family in families:
        # Decoy families (families.yaml: decoy_from_negatives) ship part of
        # their negative pool into the searchable DIAMOND db as decoy refs
        # (08a_build_decoy_refs.py, from negative.train.faa) -- simulating
        # benchmark negative reads from that same pool would let a read
        # trivially "recognize" its own literal source sequence sitting in
        # the db, inflating precision rather than testing genuine
        # discrimination. Use the disjoint negative.test.faa split instead.
        # Non-decoy families are unaffected (their negatives are never in
        # the db), so they keep using the full negative.faa as before.
        neg_fname = (f"{family}.negative.test.faa" if family in decoy_families
                     else f"{family}.negative.faa")
        for label, fname in (("positive", f"{family}.positive.test.faa"),
                              ("negative", neg_fname)):
            records = parse_fasta(args.refs / fname)
            if not records:
                continue

            rng.shuffle(records)
            held_out = records  # positive.test.faa is already held-out;
                                 # negatives use all hard negatives (or, for
                                 # decoy families, the held-out negative.test
                                 # split disjoint from what became decoys)

            genome_fasta = genome_dir / f"{family}.{label}.fasta"
            with open(genome_fasta, "w") as fh:
                for header, protein in held_out:
                    context = build_genomic_context(protein, rng)
                    fh.write(f">{header}\n{context}\n")
                    truth_rows.append({
                        "source_header": header,
                        "family": family,
                        "label": label,
                    })

            n_reads = len(held_out) * args.reads_per_sequence

            if read_lengths:
                # Exact length control needed for stratification -- wgsim
                # respects -1/-2 exactly, iss's named models don't let you
                # dial in an arbitrary bp length. Fresh seed per length (not
                # just per family/label) so multiple lengths of the same
                # sample don't draw identical read starting positions.
                for length in read_lengths:
                    out_prefix = args.out / f"{family}.{label}.rl{length}"
                    sim_seed = rng.randint(0, 2**31 - 1)
                    print(f"[{family}/{label}] {len(held_out)} sequences -> "
                          f"~{n_reads} simulated read pairs (wgsim, {length}bp, seed={sim_seed})")
                    simulate_with_wgsim(genome_fasta, out_prefix, n_reads, sim_seed, read_length=length)
            else:
                out_prefix = args.out / f"{family}.{label}"
                # Deterministic per-(family, label) seed, drawn from the
                # module's own seeded rng -- reproducible given the same
                # RANDOM_SEED and families.yaml iteration order, without
                # wgsim/iss's own internal RNGs (which neither tool seeds
                # by default) silently making every run's actual simulated
                # reads different even when the held-out input sequences
                # are identical.
                sim_seed = rng.randint(0, 2**31 - 1)
                print(f"[{family}/{label}] {len(held_out)} sequences -> "
                      f"~{n_reads} simulated read pairs ({simulator}, seed={sim_seed})")
                if simulator == "iss":
                    simulate_with_iss(genome_fasta, out_prefix, n_reads, sim_seed)
                else:
                    simulate_with_wgsim(genome_fasta, out_prefix, n_reads, sim_seed)

    if args.background:
        bg_r1, bg_r2 = args.background
        # Named/laid out exactly like every other sample ("background_R1/R2")
        # so 10_run_benchmark.sh's existing glob picks it up automatically --
        # no separate wiring needed to actually run DIAMOND/HMM against it.
        ext = ".fastq.gz" if str(bg_r1).endswith(".gz") else ".fastq"
        shutil.copy(bg_r1, args.out / f"background_R1{ext}")
        shutil.copy(bg_r2, args.out / f"background_R2{ext}")
        print(f"Copied background read pair into {args.out} as sample 'background' "
              f"(scored in background_fpr.tsv, not per-family metrics)")

    truth_path = args.out / "truth.tsv"
    with open(truth_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source_header", "family", "label"], delimiter="\t")
        writer.writeheader()
        writer.writerows(truth_rows)

    print(f"\nDone. Truth table written to {truth_path}")
    print("Next: bash 10_run_benchmark.sh")


if __name__ == "__main__":
    main()
