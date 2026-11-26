#!/usr/bin/env python3
"""
real_read_truth.py — per-read ground truth for reads simulated by wgsim
directly from a real genome, via wgsim's own embedded fragment coordinates
(no k-mer matching needed, unlike osmo_refdb's compute_read_truth.py
which has to support simulators that don't expose coordinates).

wgsim read header format (from wgsim.c):
  @<contig>_<frag_start>_<frag_end>_<e1>:<s1>:<i1>_<e2>:<s2>:<i2>_<pairhex>/<mate>
frag_start/frag_end are the 1-based outer coordinates of the simulated
fragment, the same in both mates' headers. One mate is read forward from
frag_start, the other as the reverse complement ending at frag_end, and
which mate takes which end is drawn per pair: neither the mate number nor the
name says it. Each mate is therefore placed by comparing its sequence with
the genome at both ends (forward at frag_start, reverse complement at
frag_end) and taking the end with fewer mismatches.

A read belongs to a family when it overlaps the family's gene region by at
least --min-overlap bp.

Usage:
  python real_read_truth.py --r1 reads_R1.fastq --r2 reads_R2.fastq \\
      --genome genome.fna --regions ecoli_gene_regions.tsv \\
      --out real_read_truth.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

HEADER_RE = re.compile(
    r"^(?P<contig>.+)_(?P<frag_start>\d+)_(?P<frag_end>\d+)_"
    r"\d+:\d+:\d+_\d+:\d+:\d+_[0-9a-f]+/(?P<mate>[12])$"
)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def load_regions(path: Path) -> list[tuple[str, str, int, int]]:
    regions = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            regions.append((row["family"], row["contig"], int(row["start"]), int(row["end"])))
    return regions


def load_genome(path: Path) -> dict[str, str]:
    """{contig: sequence}; contig is the header up to the first whitespace, as wgsim names it."""
    opener = gzip.open if path.suffix == ".gz" else open
    seqs, name, chunks = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def mismatches(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def read_start(seq: str, contig_seq: str, frag_start: int, frag_end: int) -> int:
    """1-based start of the read on the genome: frag_start when it reads forward
    from the fragment's left end, frag_end - len + 1 when it is the reverse
    complement of its right end (ties go to the left end)."""
    n = len(seq)
    left = contig_seq[frag_start - 1:frag_start - 1 + n]
    right = contig_seq[frag_end - n:frag_end].translate(COMPLEMENT)[::-1]
    return frag_start if mismatches(seq.upper(), left) <= mismatches(seq.upper(), right) else frag_end - n + 1


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


def iter_fastq(path: Path):
    """(read_id, sequence) per record."""
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            yield header[1:].rstrip("\n").split()[0], seq


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True, type=Path)
    ap.add_argument("--r2", required=True, type=Path)
    ap.add_argument("--genome", required=True, type=Path,
                    help="the FASTA wgsim simulated from (places each mate at its fragment end)")
    ap.add_argument("--regions", required=True, type=Path)
    ap.add_argument("--min-overlap", type=int, default=20)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    regions = load_regions(args.regions)
    genome = load_genome(args.genome)
    print(f"Loaded {len(regions)} target gene regions, {len(genome)} contigs")

    n_total = n_unparsed = n_hit = 0
    rows = []

    for fastq_path in (args.r1, args.r2):
        for read_id, seq in iter_fastq(fastq_path):
            n_total += 1
            m = HEADER_RE.match(read_id)
            if not m or m.group("contig") not in genome:
                n_unparsed += 1
                continue
            contig = m.group("contig")
            start = read_start(seq, genome[contig], int(m.group("frag_start")), int(m.group("frag_end")))
            end = start + len(seq) - 1

            family = overlapping_family(contig, start, end, regions, args.min_overlap)
            if family:
                n_hit += 1
            rows.append({
                "read_id": read_id,
                "family": family if family else "background",
            })

    if n_unparsed:
        print(f"WARNING: {n_unparsed}/{n_total} read headers didn't match the expected wgsim "
              f"format or name a contig missing from --genome -- check --r1/--r2 were produced "
              f"by wgsim from --genome")

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["read_id", "family"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{n_total} reads: {n_hit} overlap a target gene region, "
          f"{n_total - n_hit - n_unparsed} are background")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
