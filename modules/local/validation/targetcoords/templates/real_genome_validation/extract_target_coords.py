#!/usr/bin/env python3
"""
extract_target_coords.py — pull genomic coordinates for the E. coli genes
already confirmed present via `osmotool annotate` (see
ecoli_k12_both.hmmscan.tblout), from Prodigal's own coordinate-bearing
FASTA headers, so real short reads can later be checked for overlap
against these exact loci.

Usage:
  python extract_target_coords.py --prodigal-faa ecoli_k12.prodigal.faa \\
      --targets assets/real_genome_validation/target_orfs.tsv --genome ecoli_k12 \\
      --out gene_regions.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

def load_targets(path: Path, genome: str) -> dict[str, str]:
    """family -> Prodigal ORF ID for one genome, from assets/real_genome_validation/target_orfs.tsv
    (curated once against independent annotation: RefSeq gene names, UniProt)."""
    with open(path) as fh:
        return {r["family"]: r["orf_id"] for r in csv.DictReader(fh, delimiter="\t")
                if r["genome"] == genome}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prodigal-faa", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path, help="target_orfs.tsv")
    ap.add_argument("--genome", required=True, help="genome name in target_orfs.tsv")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    target_orfs = load_targets(args.targets, args.genome)
    orf_to_family = {v: k for k, v in target_orfs.items()}
    found: dict[str, tuple[str, int, int]] = {}

    for line in args.prodigal_faa.read_text().splitlines():
        if not line.startswith(">"):
            continue
        # Prodigal header: >seqid_orfnum # start # end # strand # attrs
        parts = [p.strip() for p in line[1:].split("#")]
        orf_id = parts[0].split()[0]
        if orf_id not in orf_to_family:
            continue
        family = orf_to_family[orf_id]
        start, end = int(parts[1]), int(parts[2])
        contig = orf_id.rsplit("_", 1)[0]
        found[family] = (contig, start, end)

    missing = set(target_orfs) - set(found)
    if missing:
        print(f"WARNING: couldn't find coordinates for: {sorted(missing)} "
              f"-- check --prodigal-faa is the retained file from the same "
              f"annotate run that produced ecoli_k12_both.hmmscan.tblout")

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["family", "contig", "start", "end"])
        for family, (contig, start, end) in sorted(found.items()):
            writer.writerow([family, contig, start, end])
            print(f"[{family}] {contig}:{start}-{end} ({end - start + 1}bp)")

    print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
