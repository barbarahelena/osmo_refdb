#!/usr/bin/env python3
"""
simulate_reads.py -- paired reads from real catalog genes, in nucleotides,
with each read's source gene as truth.

Each pair starts at a uniformly drawn base of the slice (so a gene gets
reads in proportion to its length), on a random strand. The sample's strain
differs from its catalog gene by --divergence per base: substitutions, and
1-3 bp indels for --indel-fraction of the changes, drawn once per gene
(seeded by the gene's index), so both mates and every read of a gene see the
same strain. Sequencing errors are substitutions at --error-rate per base.

Genes shorter than one read (150 bp) get no reads. The insert is
N(300, 30), capped at the gene's length.

Outputs: <prefix>_R1.fastq, <prefix>_R2.fastq (reads r<k>/1, r<k>/2) and
<prefix>.truth.tsv (read k -> gene, start, insert, strand).

Usage:
  python simulate_reads.py --genes slice.fna --pairs 2000000 --divergence 0.02 \\
      --error-rate 0.005 --seed 1 --out reads/d0.02
"""
from __future__ import annotations

import argparse
import bisect
import math
import random
from itertools import accumulate
from pathlib import Path

READ_LEN = 150
INSERT_MEAN, INSERT_SD = 300, 30
BASES = "ACGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def load_genes(path: Path) -> tuple[list[str], list[str]]:
    ids, seqs, seq = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq:
                    seqs.append("".join(seq).upper())
                ids.append(line[1:].split()[0])
                seq = []
            else:
                seq.append(line)
    if ids:
        seqs.append("".join(seq).upper())
    return ids, seqs


def next_event(i: int, rate: float, rng: random.Random) -> int:
    """Position of the next event at or after i in a per-base Bernoulli(rate)
    process (geometric jump), so the cost is per event, not per base."""
    return i + int(math.log(1.0 - rng.random()) / math.log(1.0 - rate))


def other_base(b: str, rng: random.Random) -> str:
    return rng.choice([x for x in BASES if x != b] or BASES)


def diverge(seq: str, rate: float, indel_fraction: float, rng: random.Random) -> str:
    """The strain's copy of a gene: about rate * len changes."""
    if rate <= 0:
        return seq
    out, i = [], 0
    while True:
        j = next_event(i, rate, rng)
        if j >= len(seq):
            out.append(seq[i:])
            return "".join(out)
        out.append(seq[i:j])
        if rng.random() < indel_fraction:
            n = rng.randint(1, 3)
            if rng.random() < 0.5:
                out.append("".join(rng.choice(BASES) for _ in range(n)) + seq[j])   # insertion
                i = j + 1
            else:
                i = j + n                                                           # deletion
        else:
            out.append(other_base(seq[j], rng))
            i = j + 1


def add_errors(read: str, rate: float, rng: random.Random) -> str:
    if rate <= 0:
        return read
    r = list(read)
    j = next_event(0, rate, rng)
    while j < len(r):
        r[j] = other_base(r[j], rng)
        j = next_event(j + 1, rate, rng)
    return "".join(r)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=Path, required=True)
    ap.add_argument("--pairs", type=int, required=True)
    ap.add_argument("--divergence", type=float, default=0.0)
    ap.add_argument("--indel-fraction", type=float, default=0.1)
    ap.add_argument("--error-rate", type=float, default=0.005)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True, help="output prefix")
    args = ap.parse_args()

    ids, seqs = load_genes(args.genes)
    weights = [max(len(s) - READ_LEN + 1, 0) for s in seqs]
    cum = list(accumulate(weights))
    total = cum[-1]
    rng = random.Random(args.seed)

    # draw every pair's gene first, so each gene's strain is made once
    picks: dict[int, list[int]] = {}
    for k in range(args.pairs):
        g = bisect.bisect_right(cum, rng.randrange(total))
        picks.setdefault(g, []).append(k)
    print(f"{len(ids):,} genes, {total:,} start positions; {args.pairs:,} pairs from {len(picks):,} genes")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    qual = "I" * READ_LEN
    pairs: list[tuple[int, str, str]] = []
    with open(f"{args.out}.truth.tsv", "w") as truth:
        truth.write("read\tgene\tstart\tinsert\tstrand\n")
        for g in sorted(picks):
            grng = random.Random(f"{args.seed}|{g}")
            strain = diverge(seqs[g], args.divergence, args.indel_fraction, grng)
            if len(strain) < READ_LEN:
                strain = seqs[g]      # a deletion-heavy strain shorter than a read: keep the gene
            for k in picks[g]:
                insert = min(max(int(rng.gauss(INSERT_MEAN, INSERT_SD)), READ_LEN), len(strain))
                start = rng.randrange(len(strain) - insert + 1)
                frag = strain[start:start + insert]
                strand = "+"
                if rng.random() < 0.5:
                    frag, strand = frag.translate(COMP)[::-1], "-"
                r1 = add_errors(frag[:READ_LEN], args.error_rate, rng)
                r2 = add_errors(frag[-READ_LEN:].translate(COMP)[::-1], args.error_rate, rng)
                pairs.append((k, r1, r2))
                truth.write(f"r{k}\t{ids[g]}\t{start + 1}\t{insert}\t{strand}\n")

    pairs.sort()   # read order independent of gene order
    with open(f"{args.out}_R1.fastq", "w") as f1, open(f"{args.out}_R2.fastq", "w") as f2:
        for k, r1, r2 in pairs:
            f1.write(f"@r{k}/1\n{r1}\n+\n{qual}\n")
            f2.write(f"@r{k}/2\n{r2}\n+\n{qual}\n")
    print(f"written: {args.out}_R1.fastq, {args.out}_R2.fastq, {args.out}.truth.tsv")


if __name__ == "__main__":
    main()
