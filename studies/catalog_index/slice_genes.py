#!/usr/bin/env python3
"""
slice_genes.py -- a seeded random slice of a gene catalog: each gene is kept
with probability n / total, so the slice holds about n genes in the
catalog's own order and proportions (about 5% with a panel HMM hit in the
human-gut catalog).

Usage:
  python slice_genes.py --genes GMGC10...fna.gz --n 1000000 --seed 1 --out slice.fna
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path


def iter_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header, seq = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=Path, required=True)
    ap.add_argument("--n", type=int, required=True, help="genes wanted (about)")
    ap.add_argument("--total", type=int, default=None,
                    help="genes in the catalog; counted in a first pass when not given")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    total = args.total or sum(1 for _ in iter_fasta(args.genes))
    p = min(1.0, args.n / total)
    rng = random.Random(args.seed)
    kept = bases = 0
    with open(args.out, "w") as out:
        for gid, seq in iter_fasta(args.genes):
            if rng.random() < p:
                out.write(f">{gid}\n{seq}\n")
                kept += 1
                bases += len(seq)
    print(f"catalog: {total:,} genes; slice: {kept:,} genes, {bases:,} bp (p = {p:.5f}) -> {args.out}")


if __name__ == "__main__":
    main()
