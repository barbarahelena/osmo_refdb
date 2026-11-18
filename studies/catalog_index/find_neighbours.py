#!/usr/bin/env python3
"""
find_neighbours.py -- the catalog genes without an HMM hit that a read could
mistake for a hit gene: the competitors a hits-only index leaves out.

  windows  every gene without a hit, cut into read-sized windows
           (--window bp every --step bp), as FASTA, for bwa-mem3 against an
           index of the hit genes
  select   from that SAM (stdin), every gene with a window aligned at
           >= --min-identity over >= --min-aligned bp (profile settles a
           read at 0.97 over 45 bp; the margin is for strain divergence and
           sequencing errors); writes the neighbours with their best identity,
           and the FASTA of hit genes plus neighbours

Identity is 1 - NM / alignment columns (M, I and D of the CIGAR), as
profile computes it. Every placement counts: the primary and those in XA.

Usage:
  python find_neighbours.py windows --genes slice.fna --hits hit_ids.txt > windows.fa
  bwa-mem3 mem -t 32 hits_idx windows.fa | \\
      python find_neighbours.py select --genes slice.fna --hits hit_ids.txt \\
          --out-tsv neighbours.tsv --out-fna hits_neighbours.fna
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def iter_fasta(path: Path):
    header, seq = None, []
    with open(path) as fh:
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


def columns_and_aligned(cigar: str) -> tuple[int, int]:
    """(alignment columns M+I+D, aligned read bases M+I)."""
    cols = aligned = 0
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in "M=X":
            cols += n
            aligned += n
        elif op == "I":
            cols += n
            aligned += n
        elif op == "D":
            cols += n
    return cols, aligned


def windows(args) -> None:
    hits = set(Path(args.hits).read_text().split())
    out = sys.stdout
    n = 0
    for gid, seq in iter_fasta(args.genes):
        if gid in hits:
            continue
        last = max(len(seq) - args.window, 0)
        starts = list(range(0, last + 1, args.step))
        if starts[-1] != last:
            starts.append(last)
        for s in starts:
            out.write(f">{gid}|{s}\n{seq[s:s + args.window]}\n")
            n += 1
    print(f"{n:,} windows", file=sys.stderr)


def select(args) -> None:
    best: dict[str, float] = {}
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & 0x4:
            continue
        gene = f[0].rsplit("|", 1)[0]
        tags = {t[:2]: t[5:] for t in f[11:]}
        placements = [(f[5], int(tags.get("NM", "0")))]
        # XA:Z:chr,pos,CIGAR,NM;...
        for alt in filter(None, tags.get("XA", "").split(";")):
            parts = alt.split(",")
            if len(parts) >= 4:
                placements.append((parts[2], int(parts[3])))
        for cigar, nm in placements:
            cols, aligned = columns_and_aligned(cigar)
            if not cols or aligned < args.min_aligned:
                continue
            ident = 1 - nm / cols
            if ident >= args.min_identity and ident > best.get(gene, 0.0):
                best[gene] = ident
    hits = set(Path(args.hits).read_text().split())
    with open(args.out_tsv, "w") as out:
        out.write("gene\tbest_identity\n")
        for g in sorted(best):
            out.write(f"{g}\t{best[g]:.4f}\n")
    keep = hits | set(best)
    n = bases = 0
    with open(args.out_fna, "w") as out:
        for gid, seq in iter_fasta(args.genes):
            if gid in keep:
                out.write(f">{gid}\n{seq}\n")
                n += 1
                bases += len(seq)
    for t in (0.90, 0.95, 0.97):
        print(f"neighbours at >= {t:.2f}: {sum(1 for v in best.values() if v >= t):,}", file=sys.stderr)
    print(f"{len(hits):,} hit genes + {len(best):,} neighbours = {n:,} genes, {bases:,} bp -> {args.out_fna}",
          file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("windows")
    w.add_argument("--genes", type=Path, required=True)
    w.add_argument("--hits", type=Path, required=True, help="IDs of the genes with an HMM hit, one per line")
    w.add_argument("--window", type=int, default=150)
    w.add_argument("--step", type=int, default=50)
    s = sub.add_parser("select")
    s.add_argument("--genes", type=Path, required=True)
    s.add_argument("--hits", type=Path, required=True)
    s.add_argument("--min-identity", type=float, default=0.90)
    s.add_argument("--min-aligned", type=int, default=45)
    s.add_argument("--out-tsv", type=Path, required=True)
    s.add_argument("--out-fna", type=Path, required=True)
    args = ap.parse_args()
    windows(args) if args.cmd == "windows" else select(args)


if __name__ == "__main__":
    main()
