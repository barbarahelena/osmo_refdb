#!/usr/bin/env python3
"""
real_read_truth.py — per-read ground truth for reads simulated by wgsim
directly from a real genome, via wgsim's own embedded read coordinates
(no k-mer matching needed, unlike osmo_refdb's 09b_compute_read_truth.py
which has to support simulators that don't expose coordinates).

wgsim read header format (from wgsim.c):
  @<contig>_<start1>_<start2>_<e1>:<s1>:<i1>_<e2>:<s2>:<i2>_<pairhex>/<mate>
start1/start2 are 1-based genomic start coordinates for mate 1 and mate 2
respectively; both mates' headers carry both coordinates, differing only
in the trailing /1 or /2.

Usage:
  python real_read_truth.py --r1 reads_R1.fastq --r2 reads_R2.fastq \\
      --regions ecoli_gene_regions.tsv --read-length 150 \\
      --out real_read_truth.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

HEADER_RE = re.compile(
    r"^(?P<contig>.+)_(?P<start1>\d+)_(?P<start2>\d+)_"
    r"\d+:\d+:\d+_\d+:\d+:\d+_[0-9a-f]+/(?P<mate>[12])$"
)


def load_regions(path: Path) -> list[tuple[str, str, int, int]]:
    regions = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            regions.append((row["family"], row["contig"], int(row["start"]), int(row["end"])))
    return regions


def overlapping_family(
    contig: str, read_start: int, read_end: int,
    regions: list[tuple[str, str, int, int]], min_overlap: int,
) -> str | None:
    for family, region_contig, region_start, region_end in regions:
        if contig != region_contig:
            continue
        overlap = min(read_end, region_end) - max(read_start, region_start) + 1
        if overlap >= min_overlap:
            return family
    return None


def iter_fastq_ids(path: Path):
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            fh.readline()
            fh.readline()
            fh.readline()
            yield header[1:].rstrip("\n").split()[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True, type=Path)
    ap.add_argument("--r2", required=True, type=Path)
    ap.add_argument("--regions", required=True, type=Path)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--min-overlap", type=int, default=20)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    regions = load_regions(args.regions)
    print(f"Loaded {len(regions)} target gene regions")

    n_total = n_unparsed = n_hit = 0
    rows = []

    for fastq_path in (args.r1, args.r2):
        for read_id in iter_fastq_ids(fastq_path):
            n_total += 1
            m = HEADER_RE.match(read_id)
            if not m:
                n_unparsed += 1
                continue
            contig = m.group("contig")
            mate = m.group("mate")
            start = int(m.group("start1")) if mate == "1" else int(m.group("start2"))
            end = start + args.read_length - 1

            family = overlapping_family(contig, start, end, regions, args.min_overlap)
            if family:
                n_hit += 1
            rows.append({
                "read_id": read_id,
                "family": family if family else "background",
            })

    if n_unparsed:
        print(f"WARNING: {n_unparsed}/{n_total} read headers didn't match the "
              f"expected wgsim format -- check --r1/--r2 were actually produced by wgsim")

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["read_id", "family"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{n_total} reads: {n_hit} overlap a target gene region, "
          f"{n_total - n_hit - n_unparsed} are background")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
