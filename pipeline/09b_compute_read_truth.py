#!/usr/bin/env python3
"""
09b_compute_read_truth.py — per-read ground truth via actual CDS overlap.

Motivating bug: 09_simulate_reads.py's truth.tsv labels every read simulated
from a positive/negative construct with that construct's label, but each
construct is FLANK_LEN (300bp) of random DNA + the real CDS + FLANK_LEN more
random DNA. A read landing entirely in the flanking DNA carries no gene
signal at all -- scoring it as a "positive" read that a caller must recall
is wrong; a caller that correctly ignores it is being penalized as if it
missed a true positive.

This step locates each simulated read on its source contig via a k-mer seed
match (works regardless of which simulator produced the reads -- neither
InSilicoSeq nor wgsim reliably exposes read coordinates through the CLI used
in 09_simulate_reads.py) and checks overlap against the known CDS interval
[FLANK_LEN, contig_len - FLANK_LEN). Reads that can't be confidently located
(should be rare -- only at very high simulated error rates) are dropped
rather than guessed at.

Also picks up real background reads (09_simulate_reads.py --background,
written as background_R1/R2.fastq) -- these have no genomic-context contig
to check overlap against, so every one is recorded as a trivial true
negative rather than k-mer-matched. And records each read's length, so
11_compute_metrics.py can report metrics stratified by read length even
when 09_simulate_reads.py --read-lengths simulated several.

Usage:
  python 09b_compute_read_truth.py --reads results/reads \\
      --min-cds-overlap 20
Output:
  results/reads/read_truth.tsv  (read_id, family, label, uniprot_id,
                                  overlaps_cds, overlap_bp, read_length)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

FLANK_LEN = 300  # must match 09_simulate_reads.py:FLANK_LEN
KMER_SIZE = 21
N_PROBES = 8
HEADER_ID_RE = re.compile(r"\|([A-Z0-9]+)\|")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def iter_fastq_ids_and_seqs(path: Path):
    with open_maybe_gz(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()  # '+'
            fh.readline()  # quality
            read_id = header[1:].rstrip("\n").split()[0]
            yield read_id, seq


def build_kmer_index(seq: str, k: int = KMER_SIZE) -> dict[str, int]:
    return {seq[i:i + k]: i for i in range(len(seq) - k + 1)}


def locate_read(read_seq: str, contig_seq: str, kmer_index: dict[str, int],
                 k: int = KMER_SIZE, n_probes: int = N_PROBES) -> tuple[int, int] | None:
    """Return the read's (start, end) 0-based half-open interval on contig_seq,
    or None if no confident k-mer anchor was found in either orientation."""
    for candidate in (read_seq, revcomp(read_seq)):
        L = len(candidate)
        if L < k:
            continue
        step = max(1, (L - k) // n_probes)
        for i in range(0, L - k + 1, step):
            pos = kmer_index.get(candidate[i:i + k])
            if pos is None:
                continue
            start, end = pos - i, pos - i + L
            if 0 <= start and end <= len(contig_seq):
                return start, end
    return None


def find_source_contig(read_id: str, contig_headers: list[str]) -> str | None:
    """Read IDs are the contig header + simulator-appended suffix (e.g.
    '<header>_0_0/1'), so the source contig is whichever known header the
    read ID starts with. Pick the longest match in case one header happens
    to be a prefix of another."""
    matches = [h for h in contig_headers if read_id.startswith(h)]
    if not matches:
        return None
    return max(matches, key=len)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", type=Path, default=Path("results/reads"))
    ap.add_argument("--min-cds-overlap", type=int, default=20,
                     help="Minimum bp of read/CDS overlap to count a read as "
                          "carrying real gene signal")
    args = ap.parse_args()

    genome_dir = args.reads / "genomic_context"
    out_rows = []
    n_total = n_unresolved = n_overlapping = 0

    for genome_fasta in sorted(genome_dir.glob("*.fasta")):
        stem = genome_fasta.stem  # "<family>.<label>"
        family, _, label = stem.partition(".")
        contigs = parse_fasta(genome_fasta)
        contig_headers = list(contigs.keys())
        cds_bounds = {}
        for header, seq in contigs.items():
            cds_start, cds_end = FLANK_LEN, len(seq) - FLANK_LEN
            if cds_end <= cds_start:
                print(f"[{stem}] WARNING: contig '{header}' shorter than "
                      f"2*FLANK_LEN, skipping CDS bounds")
                continue
            cds_bounds[header] = (cds_start, cds_end)
        kmer_indexes: dict[str, dict[str, int]] = {}

        # Glob rather than a fixed filename: read-length-stratified runs
        # (09_simulate_reads.py --read-lengths) write e.g.
        # "<family>.<label>.rl150_R1.fastq" instead of "<family>.<label>_R1.fastq".
        for r1_path in sorted(args.reads.glob(f"{stem}*_R1.fastq*")):
            r2_path = Path(str(r1_path).replace("_R1.fastq", "_R2.fastq"))
            if not r2_path.exists():
                continue

            for fastq_path in (r1_path, r2_path):
                for read_id, read_seq in iter_fastq_ids_and_seqs(fastq_path):
                    n_total += 1
                    contig_header = find_source_contig(read_id, contig_headers)
                    if contig_header is None or contig_header not in cds_bounds:
                        n_unresolved += 1
                        continue

                    contig_seq = contigs[contig_header]
                    if contig_header not in kmer_indexes:
                        kmer_indexes[contig_header] = build_kmer_index(contig_seq)

                    located = locate_read(read_seq, contig_seq, kmer_indexes[contig_header])
                    if located is None:
                        n_unresolved += 1
                        continue

                    read_start, read_end = located
                    cds_start, cds_end = cds_bounds[contig_header]
                    overlap_bp = max(0, min(read_end, cds_end) - max(read_start, cds_start))
                    overlaps_cds = overlap_bp >= args.min_cds_overlap
                    n_overlapping += int(overlaps_cds)

                    m = HEADER_ID_RE.search(contig_header)
                    uniprot_id = m.group(1) if m else contig_header

                    out_rows.append({
                        "read_id": read_id,
                        "family": family,
                        "label": label,
                        "uniprot_id": uniprot_id,
                        "overlaps_cds": int(overlaps_cds),
                        "overlap_bp": overlap_bp,
                        "read_length": len(read_seq),
                    })

    # Background reads (09_simulate_reads.py --background): real reads with
    # no known osmoadaptation content and no genomic-context contig to check
    # overlap against, so every one is trivially a true negative for every
    # family -- no k-mer matching needed.
    for r1_path in sorted(args.reads.glob("background_R1.fastq*")):
        r2_path = Path(str(r1_path).replace("_R1.fastq", "_R2.fastq"))
        if not r2_path.exists():
            continue
        for fastq_path in (r1_path, r2_path):
            for read_id, read_seq in iter_fastq_ids_and_seqs(fastq_path):
                n_total += 1
                out_rows.append({
                    "read_id": read_id,
                    "family": "background",
                    "label": "negative",
                    "uniprot_id": "",
                    "overlaps_cds": 0,
                    "overlap_bp": 0,
                    "read_length": len(read_seq),
                })

    if not out_rows:
        print("No reads processed -- check --reads points at a populated "
              "results/reads directory (run 09_simulate_reads.py first).")
        return

    out_path = args.reads / "read_truth.tsv"
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    n_resolved = n_total - n_unresolved
    print(f"\n{n_total} reads seen, {n_resolved} located on their source contig "
          f"({n_unresolved} unresolved and dropped).")
    print(f"Of located reads: {n_overlapping} overlap the true CDS by "
          f">= {args.min_cds_overlap}bp; {n_resolved - n_overlapping} are "
          f"flanking-DNA-only reads with no real gene signal.")
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()