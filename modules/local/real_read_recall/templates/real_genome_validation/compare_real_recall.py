#!/usr/bin/env python3
"""
compare_real_recall.py — DIAMOND vs HMM precision/recall on reads
simulated directly from a real genome, against coordinate-based ground
truth (real_read_truth.py), instead of osmo_refdb's synthetic
UniProt-derived benchmark reads.

Usage:
  python compare_real_recall.py --truth real_read_truth.tsv \\
      --diamond-blastx ecoli_reads.blastx.tsv \\
      --hmm-tblout ecoli_reads.hmmscan.tblout \\
      --out real_recall_comparison.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ORF_SUFFIX_RE = re.compile(r"^(.*/[12])_\d+_\d+_\d+$")


def load_truth(path: Path) -> dict[str, str]:
    truth = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = row["family"]
    return truth


def load_diamond_best_hits(path: Path) -> dict[str, tuple[str, float]]:
    """{read_id: (family, bitscore)} keeping only the best hit per read."""
    best: dict[str, tuple[str, float]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            qseqid, sseqid, bitscore = parts[0], parts[1], float(parts[-1])
            family = sseqid.split("|")[0]
            if qseqid not in best or bitscore > best[qseqid][1]:
                best[qseqid] = (family, bitscore)
    return best


def load_hmm_best_hits(path: Path) -> dict[str, tuple[str, float]]:
    """{read_id: (family, bitscore)}, ORF id mapped back to read id."""
    best: dict[str, tuple[str, float]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            family, orf_id, bitscore = fields[0], fields[2], float(fields[5])
            m = ORF_SUFFIX_RE.match(orf_id)
            read_id = m.group(1) if m else orf_id
            if read_id not in best or bitscore > best[read_id][1]:
                best[read_id] = (family, bitscore)
    return best


def summarize(truth: dict[str, str], calls: dict[str, tuple[str, float]],
              families: list[str]) -> list[dict]:
    rows = []
    for family in families:
        tp = fp = fn = 0
        for read_id, true_family in truth.items():
            called = calls.get(read_id, (None, None))[0]
            is_true_positive_source = (true_family == family)
            if called == family and is_true_positive_source:
                tp += 1
            elif called == family and not is_true_positive_source:
                fp += 1
            elif called != family and is_true_positive_source:
                fn += 1
        precision = tp / (tp + fp) if (tp + fp) else float("nan")
        recall = tp / (tp + fn) if (tp + fn) else float("nan")
        rows.append({
            "family": family, "tp": tp, "fp": fp, "fn": fn,
            "precision": round(precision, 3) if precision == precision else "",
            "recall": round(recall, 3) if recall == recall else "",
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--diamond-blastx", required=True, type=Path)
    ap.add_argument("--hmm-tblout", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    truth = load_truth(args.truth)
    families = sorted({f for f in truth.values() if f != "background"})
    print(f"Ground truth: {len(truth)} reads, target families: {families}")

    diamond_hits = load_diamond_best_hits(args.diamond_blastx)
    hmm_hits = load_hmm_best_hits(args.hmm_tblout)
    print(f"DIAMOND: {len(diamond_hits)} reads with >=1 hit")
    print(f"HMM (raw score): {len(hmm_hits)} reads with >=1 hit")

    diamond_rows = summarize(truth, diamond_hits, families)
    hmm_rows = summarize(truth, hmm_hits, families)

    all_rows = [{"method": "diamond", **r} for r in diamond_rows] + \
               [{"method": "hmm", **r} for r in hmm_rows]

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n{'family':<8} {'diamond_recall':>15} {'hmm_recall':>12} {'diamond_prec':>13} {'hmm_prec':>10}")
    hmm_by_fam = {r["family"]: r for r in hmm_rows}
    for r in diamond_rows:
        h = hmm_by_fam[r["family"]]
        print(f"{r['family']:<8} {str(r['recall']):>15} {str(h['recall']):>12} "
              f"{str(r['precision']):>13} {str(h['precision']):>10}")

    print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
