#!/usr/bin/env python3
"""
03_split_train_test.py — split each family's positive FASTA into train/test
sets *before* building the HMM or DIAMOND database, so the benchmark's
held-out set is genuinely unseen by both methods (no data leakage).

Reads refs/<family>.positive.faa and writes:
  refs/<family>.positive.train.faa   (used to build HMM + DIAMOND db)
  refs/<family>.positive.test.faa    (held out, used only for benchmarking)

Also splits every family's refs/<family>.negative.faa, writing
refs/<family>.negative.train.faa / refs/<family>.negative.test.faa:

  negative.train -> cutoff calibration (05/06 HMM GA, 08b DIAMOND bitscore)
                    and, for `decoy_from_negatives` families, the decoy refs
                    built by 08a_build_decoy_refs.py
  negative.test  -> benchmark negative reads (09_simulate_reads.py), never
                    seen by calibration or the db

Before this, only decoy families split their negatives; every other family
calibrated its cutoffs on the full negative.faa and then simulated its
benchmark negative reads from that same file. `osmotool profile` applies
the 08b DIAMOND cutoffs during the benchmark, so DIAMOND precision was
measured partly on the negatives that had set the cutoff.

Decoy families keep --test-fraction (0.2): their negative.train becomes the
decoy reference set, so it stays large. Other families use
--negative-test-fraction (default 0.5): their negative.train only feeds
calibration, and the benchmark's false-positive estimate needs enough
negative reads. Families with fewer than MIN_NEG_FOR_SPLIT negatives are not
split: no split files are written (stale ones are removed), and downstream
steps fall back to the full negative.faa with a warning, i.e. the old
behaviour.

Negatives are always split at random, also under --split-mode taxonomy:
they are other genes, so holding out genera tests nothing extra there.

Two split modes:
  random (default) -- shuffle all sequences and split by test-fraction.
    Held-out sequences are still drawn from the same overall pool as
    training, so this measures generalization to unseen *sequences*, not
    unseen *lineages* -- a close relative of every test sequence is likely
    still in train.
  taxonomy -- hold out whole genera (parsed from the "tag|uniprot_id|
    organism" FASTA header this pipeline already writes) so no sequence
    from a held-out genus is seen during training. This is the scenario
    where profile HMMs classically have an edge over pairwise identity
    search (detecting something meaningfully more divergent than anything
    trained on), which a random split doesn't stress at all.

Usage:
  python 03_split_train_test.py --refs refs --test-fraction 0.2 \\
      --families families.yaml --split-mode random|taxonomy
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import yaml

RANDOM_SEED = 42
MIN_NEG_FOR_SPLIT = 5


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def load_decoy_family_names(path: Path) -> set[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return {fam["name"] for fam in data["families"] if fam.get("decoy_from_negatives")}


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


def genus_of(header: str) -> str:
    """First token of the organism name in a 'tag|uniprot_id|organism' header
    (retag_fasta's format), i.e. the genus in standard binomial nomenclature."""
    parts = header.split("|")
    organism = parts[2] if len(parts) >= 3 else ""
    return organism.split("_")[0] if organism else "unknown"


def split_random(records: list[tuple[str, str]], test_fraction: float,
                  rng: random.Random) -> tuple[list, list]:
    records = list(records)
    rng.shuffle(records)
    n_test = max(1, int(len(records) * test_fraction))
    return records[n_test:], records[:n_test]  # train, test


def split_by_taxonomy(records: list[tuple[str, str]], test_fraction: float,
                       rng: random.Random) -> tuple[list, list, int]:
    """Hold out whole genera rather than individual sequences. Returns
    (train, test, n_genera_held_out). Falls back to nothing held out (all
    train) if there's only one genus -- taxonomy-based holdout is meaningless
    with a single genus, and the caller should fall back to a random split."""
    by_genus: dict[str, list[tuple[str, str]]] = {}
    for header, seq in records:
        by_genus.setdefault(genus_of(header), []).append((header, seq))

    genera = list(by_genus.keys())
    rng.shuffle(genera)

    target_n_test = max(1, int(len(records) * test_fraction))
    test, train = [], []
    n_test_so_far = 0
    held_out_genera = 0
    for genus in genera:
        genus_records = by_genus[genus]
        # Keep at least one genus in train no matter what, so the HMM/DIAMOND
        # db is never built from zero sequences.
        if n_test_so_far < target_n_test and held_out_genera < len(genera) - 1:
            test.extend(genus_records)
            n_test_so_far += len(genus_records)
            held_out_genera += 1
        else:
            train.extend(genus_records)
    return train, test, held_out_genera


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", required=True, type=Path)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--test-fraction", type=float, default=0.2)
    ap.add_argument("--negative-test-fraction", type=float, default=0.5,
                     help="Fraction of each non-decoy family's negatives held "
                          "out for benchmark reads; the rest calibrates cutoffs. "
                          "Decoy families use --test-fraction instead.")
    ap.add_argument("--split-mode", choices=["random", "taxonomy"], default="random",
                     help="'taxonomy' holds out whole genera instead of "
                          "individual sequences, to test detection of "
                          "genuinely divergent (not just held-out) homologs")
    args = ap.parse_args()

    families = load_family_names(args.families)
    decoy_families = load_decoy_family_names(args.families)

    for family in families:
        # Independent per-family RNG (seeded from RANDOM_SEED + family name),
        # not one shared stream walked across every family in order. A
        # single shared stream means any change to how many random draws an
        # EARLIER family consumes (e.g. adding a negative split for a decoy
        # family) silently reshuffles every LATER family's train/test split
        # too -- observed in practice: enabling betL's decoy split shifted
        # otsB's and galE's benchmark numbers by real amounts even though
        # neither family's own logic changed. Per-family seeding makes each
        # family's split depend only on its own data, never on what
        # families.yaml lists before it or how other families are configured.
        rng = random.Random(f"{RANDOM_SEED}:{family}")

        pos_path = args.refs / f"{family}.positive.faa"
        records = parse_fasta(pos_path)
        if len(records) < 5:
            print(f"[{family}] SKIP: only {len(records)} positive sequences, "
                  f"too few for a train/test split")
            continue

        if args.split_mode == "taxonomy":
            n_genera = len({genus_of(h) for h, _ in records})
            if n_genera < 2:
                print(f"[{family}] only 1 genus represented -- taxonomy split "
                      f"is meaningless here, falling back to random split")
                train_records, test_records = split_random(records, args.test_fraction, rng)
            else:
                train_records, test_records, n_held_out = split_by_taxonomy(
                    records, args.test_fraction, rng)
                print(f"[{family}] held out {n_held_out}/{n_genera} genera for testing")
        else:
            train_records, test_records = split_random(records, args.test_fraction, rng)

        write_fasta(args.refs / f"{family}.positive.train.faa", train_records)
        write_fasta(args.refs / f"{family}.positive.test.faa", test_records)

        print(f"[{family}] {len(records)} total -> "
              f"{len(train_records)} train / {len(test_records)} test")

        # Negatives: every family now, not just decoy families (see module
        # docstring). Same per-family rng, drawn after the positive split, so
        # positive splits are unchanged by this and decoy families' negative
        # splits are byte-identical to before.
        neg_train_path = args.refs / f"{family}.negative.train.faa"
        neg_test_path = args.refs / f"{family}.negative.test.faa"
        neg_records = parse_fasta(args.refs / f"{family}.negative.faa")
        if len(neg_records) < MIN_NEG_FOR_SPLIT:
            for stale in (neg_train_path, neg_test_path):
                stale.unlink(missing_ok=True)
            print(f"[{family}] WARNING: only {len(neg_records)} negative sequences, "
                  f"too few to split -- calibration and benchmark will both "
                  f"use the full negative.faa")
            continue
        is_decoy = family in decoy_families
        neg_fraction = args.test_fraction if is_decoy else args.negative_test_fraction
        neg_train, neg_test = split_random(neg_records, neg_fraction, rng)
        write_fasta(neg_train_path, neg_train)
        write_fasta(neg_test_path, neg_test)
        role = "decoy refs + calibration" if is_decoy else "calibration"
        print(f"[{family}] {len(neg_records)} negatives -> "
              f"{len(neg_train)} negative.train ({role}) / "
              f"{len(neg_test)} negative.test (benchmark reads)")

    print("\nDone. Use *.positive.train.faa to build HMMs/DIAMOND db, "
          "*.positive.test.faa as held-out benchmark positives, "
          "*.negative.train.faa for calibration and *.negative.test.faa "
          "for benchmark negatives.")


if __name__ == "__main__":
    main()
