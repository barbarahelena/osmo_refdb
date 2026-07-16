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


def simulate_with_iss(genome_fasta: Path, out_prefix: Path, n_reads: int) -> None:
    subprocess.run(
        ["iss", "generate", "--genomes", str(genome_fasta),
         "--n_reads", str(n_reads), "--model", "hiseq",
         "--output", str(out_prefix)],
        check=True,
    )


def simulate_with_wgsim(genome_fasta: Path, out_prefix: Path, n_reads: int) -> None:
    # wgsim needs per-sequence coverage; simplest robust approach is to
    # concatenate all contigs and let wgsim sample reads genome-wide.
    r1 = f"{out_prefix}_R1.fastq"
    r2 = f"{out_prefix}_R2.fastq"
    n_pairs = max(1, n_reads // 2)
    subprocess.run(
        ["wgsim", "-N", str(n_pairs), "-1", "150", "-2", "150",
         str(genome_fasta), r1, r2],
        check=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs", required=True, type=Path, help="Path to osmo_refdb/refs")
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--reads-per-sequence", type=int, default=20,
                     help="Simulated read pairs per held-out sequence")
    ap.add_argument("--background", nargs="*", default=None,
                     help="Optional real background FASTQ file(s) with no "
                          "known osmoadaptation content, copied in as extra "
                          "true negatives for false-positive rate estimation")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rng = random.Random(RANDOM_SEED)
    families = load_family_names(args.families)
    args.out.mkdir(parents=True, exist_ok=True)
    genome_dir = args.out / "genomic_context"
    genome_dir.mkdir(exist_ok=True)

    simulator = which_simulator()

    truth_rows = []

    for family in families:
        for label, fname in (("positive", f"{family}.positive.test.faa"),
                              ("negative", f"{family}.negative.faa")):
            records = parse_fasta(args.refs / fname)
            if not records:
                continue

            rng.shuffle(records)
            held_out = records  # positive.test.faa is already held-out;
                                 # negatives use all hard negatives

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
            out_prefix = args.out / f"{family}.{label}"
            print(f"[{family}/{label}] {len(held_out)} sequences -> "
                  f"~{n_reads} simulated read pairs ({simulator})")
            if simulator == "iss":
                simulate_with_iss(genome_fasta, out_prefix, n_reads)
            else:
                simulate_with_wgsim(genome_fasta, out_prefix, n_reads)

    if args.background:
        bg_dir = args.out / "background"
        bg_dir.mkdir(exist_ok=True)
        for bg_file in args.background:
            shutil.copy(bg_file, bg_dir / Path(bg_file).name)
        print(f"Copied {len(args.background)} background FASTQ file(s) into {bg_dir}")

    truth_path = args.out / "truth.tsv"
    with open(truth_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source_header", "family", "label"], delimiter="\t")
        writer.writeheader()
        writer.writerows(truth_rows)

    print(f"\nDone. Truth table written to {truth_path}")
    print("Next: bash 10_run_benchmark.sh")


if __name__ == "__main__":
    main()
