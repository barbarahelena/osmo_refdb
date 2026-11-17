#!/usr/bin/env python3
"""
08a_build_decoy_refs.py — turn a family's QC'd negative set into searchable
DIAMOND decoy references, for families where a single score threshold
cannot separate the target from a specific, named confusable paralog.

Motivating case: betL vs betT/caiT (BCCT-family transporters). v4's
calibration data showed betL's positive and negative score distributions
overlap almost completely (min_positive_score=87.2 vs
max_negative_score=642.1 in hmms/cutoff_manifest.tsv) -- no absolute
DIAMOND bitscore or HMM GA cutoff can cleanly separate them. betL's
negative_query already targets its only two real confusable paralogs by
name (betT, caiT), so after 01b/01c QC, refs/betL.negative.faa IS a clean,
purpose-built decoy set -- it just isn't in the searchable database yet.

Reads refs/<family>.negative.train.faa (written by 03_split_train_test.py,
only for families marked decoy_from_negatives), NOT the full
refs/<family>.negative.faa -- 09_simulate_reads.py simulates benchmark
negative reads from the disjoint refs/<family>.negative.test.faa split for
these families specifically, so a benchmark read can never trivially
"recognize" its own literal source sequence sitting in the db as a decoy.
Run 03_split_train_test.py before this script.

This step doesn't add any new scoring logic: it relies entirely on
quantifier.select_best_hits(), which already picks the single
highest-scoring alignment across the *whole* combined DIAMOND db. Today
betT/caiT sequences aren't in that db, so they can never win a read/protein
away from a weak betL hit. Adding them as "<family>_decoy"-labelled entries
makes that existing best-hit contest competitive for real betT/caiT
signal, instead of leaning on a threshold the calibration data proves can't
work. Decoy hits are excluded from all reported output (see
08c_write_scope_manifest.py) -- they exist purely to be won against, never
to be counted.

Only families explicitly marked `decoy_from_negatives: true` in
families.yaml are processed -- this only helps when negative_query targets
specific named paralogs (as betL's does); a broad Pfam pool would make a
weak, noisy decoy and shouldn't be used this way.

Dedup against positive references: a decoy candidate is dropped if its
accession is already used as a positive training reference for ANY family
in the combined DIAMOND db (not just this decoy family's own siblings).
Without this, a decoy set built from named-paralog negatives can contain
the exact same sequence that's simultaneously a real positive reference
elsewhere (e.g. a proX_decoy entry that's byte-identical to opuAC's own
positive reference) -- select_best_hits() could then tie between the
correct positive call and the decoy-labelled duplicate, and since decoys
are excluded from reported output, an unlucky tie silently drops the read
instead of correctly crediting the real family. Filtering these out by
construction closes the collision at the source rather than relying on a
benchmark to catch it after the fact.

Usage:
  python 08a_build_decoy_refs.py --refs refs --families families.yaml
Output:
  refs/<family>.decoy.faa   (only for families with decoy_from_negatives: true)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def load_decoy_families(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"] if fam.get("decoy_from_negatives")]


def load_all_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def header_accession(header: str) -> str:
    """'<family>|UniProtID|Organism' -> 'UniProtID'."""
    parts = header.split("|")
    return parts[1] if len(parts) > 1 else header


def load_all_positive_accessions(refs_dir: Path, family_names: list[str]) -> set[str]:
    """Accessions already used as a positive training reference for any
    family in the panel, regardless of whether that family builds decoys."""
    accessions: set[str] = set()
    for family in family_names:
        pos_path = refs_dir / f"{family}.positive.train.faa"
        for header, _ in parse_fasta(pos_path):
            accessions.add(header_accession(header))
    return accessions


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header, chunks = None, []
    if not path.exists():
        return records
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def retag_as_decoy(header: str, family: str) -> str:
    """<anything>|UniProtID|Organism -> <family>_decoy|UniProtID|Organism.
    Only the first '|'-field (gene_family_from_header's family label) is
    replaced; the rest of the header (UniProt ID, organism) is preserved
    for traceability."""
    parts = header.split("|", 1)
    rest = parts[1] if len(parts) > 1 else ""
    return f"{family}_decoy|{rest}" if rest else f"{family}_decoy"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    args = ap.parse_args()

    decoy_families = load_decoy_families(args.families)
    if not decoy_families:
        print("No families marked decoy_from_negatives -- nothing to build.")
        return

    all_family_names = load_all_family_names(args.families)
    positive_accessions = load_all_positive_accessions(args.refs, all_family_names)

    for family in decoy_families:
        neg_path = args.refs / f"{family}.negative.train.faa"
        decoy_path = args.refs / f"{family}.decoy.faa"

        if not neg_path.exists() or neg_path.stat().st_size == 0:
            print(f"[{family}] SKIP: no negative train split found at {neg_path} "
                  f"-- run 01/01b/01c then 03_split_train_test.py first")
            continue

        records = parse_fasta(neg_path)
        n_dup = sum(1 for header, _ in records if header_accession(header) in positive_accessions)
        records = [(h, s) for h, s in records if header_accession(h) not in positive_accessions]
        retagged = [(retag_as_decoy(header, family), seq) for header, seq in records]
        write_fasta(decoy_path, retagged)
        dup_note = f" ({n_dup} dropped as duplicate positive references elsewhere)" if n_dup else ""
        print(f"[{family}] {len(retagged)} decoy references written -> {decoy_path}{dup_note}")

    print("\nDone.")


if __name__ == "__main__":
    main()
