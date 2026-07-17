#!/usr/bin/env python3
"""
extract_target_coords.py — pull genomic coordinates for the E. coli genes
already confirmed present via `osmotool annotate` (see
ecoli_k12_both.hmmscan.tblout), from Prodigal's own coordinate-bearing
FASTA headers, so real short reads can later be checked for overlap
against these exact loci.

Usage:
  python extract_target_coords.py --prodigal-faa ecoli_k12.prodigal.faa \\
      --out ecoli_gene_regions.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

# family -> Prodigal ORF ID, from ecoli_k12_both.hmmscan.tblout (the HMM
# calls we validated against known E. coli biology). ectA excluded on
# purpose: it's the known false positive, not a real gene to search for.
TARGET_ORFS = {
    "nhaA": "NC_000913.3_17",
    "kdpA": "NC_000913.3_677",
    "otsA": "NC_000913.3_1873",
    "otsB": "NC_000913.3_1874",
    "proX": "NC_000913.3_2633",
    "mscS": "NC_000913.3_2871",
    "mscL": "NC_000913.3_3219",
    "proP": "NC_000913.3_4029",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prodigal-faa", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    orf_to_family = {v: k for k, v in TARGET_ORFS.items()}
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

    missing = set(TARGET_ORFS) - set(found)
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
