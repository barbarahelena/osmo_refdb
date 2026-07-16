#!/usr/bin/env python3
"""
03_split_train_test.py — split each family's positive FASTA into train/test
sets *before* building the HMM or DIAMOND database, so the benchmark's
held-out set is genuinely unseen by both methods (no data leakage).

Reads refs/<family>.positive.faa and writes:
  refs/<family>.positive.train.faa   (used to build HMM + DIAMOND db)
  refs/<family>.positive.test.faa    (held out, used only for benchmarking)

Usage:
  python 03_split_train_test.py --refs refs --test-fraction 0.2 --families families.yaml
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import yaml

RANDOM_SEED = 42


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def parse_fasta(path: Path) -> list[tuple[str, str]]:
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


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs", required=True, type=Path)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--test-fraction", type=float, default=0.2)
    args = ap.parse_args()

    rng = random.Random(RANDOM_SEED)
    families = load_family_names(args.families)

    for family in families:
        pos_path = args.refs / f"{family}.positive.faa"
        records = parse_fasta(pos_path)
        if len(records) < 5:
            print(f"[{family}] SKIP: only {len(records)} positive sequences, "
                  f"too few for a train/test split")
            continue

        rng.shuffle(records)
        n_test = max(1, int(len(records) * args.test_fraction))
        test_records = records[:n_test]
        train_records = records[n_test:]

        write_fasta(args.refs / f"{family}.positive.train.faa", train_records)
        write_fasta(args.refs / f"{family}.positive.test.faa", test_records)

        print(f"[{family}] {len(records)} total -> "
              f"{len(train_records)} train / {len(test_records)} test")

    print("\nDone. Use *.positive.train.faa to build HMMs/DIAMOND db, "
          "*.positive.test.faa as held-out benchmark positives.")


if __name__ == "__main__":
    main()
