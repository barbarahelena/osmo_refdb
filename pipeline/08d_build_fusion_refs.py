#!/usr/bin/env python3
"""
08d_build_fusion_refs.py — turn fusion-candidate sequences flagged by
01c_check_length_outliers.py into a single, searchable DIAMOND reference
per declared fusion pair (families.yaml: `fusion_partner`).

Motivating case: mrpA/mrpB (Task 1b). In some lineages the Mrp/Mnh Na+/H+
antiporter's A and B subunits occur as one fused ORF instead of two genes,
confirmed via direct UniProt search and dispersed across Actinomycetota,
Bacillota (specifically Paenibacillaceae), and some Alphaproteobacteria --
see families.yaml's mrpA/mrpB entries for the taxonomic evidence. 01c
already keeps these out of the standard positive.faa (a ~1000aa fused
sequence would badly gap the MAFFT alignment of ~800aa standalone mrpA
sequences) and routes them to refs/<family>.positive.fusion_candidates.faa
instead of the length-outliers pile. This script is what actually makes
them searchable: it merges a pair's fusion candidates (deduped by UniProt
accession, since the same fused protein can turn up in both mrpA's and
mrpB's raw fetch depending which gene symbol it happened to be annotated
under), retags them under one shared "<familyA>_<familyB>_fused" label, and
writes a combined FASTA that 08_build_diamond_db.sh folds into the DIAMOND
db like any other family.

Unlike 08a's decoy references, a fused-ORF hit is a REAL, reportable
detection target (it genuinely carries both subunits), not a sink for
mislabeled calls -- it is NOT added to the profile/annotate exclusion lists
in 08c_write_scope_manifest.py. Complex-aware scoring for what a single
"<familyA>_<familyB>_fused" hit implies for both individual family calls
lives downstream, in osmotool (Task 2), not in this repo.

Each family in a fusion_partner pair processes the pair once (the second
family sees the label already built and skips it), so this is safe however
many pairs are declared.

KNOWN LIMITATION: 08b_calibrate_diamond_cutoffs.py calibrates a per-family
minimum bitscore only for names literally listed in families.yaml -- a
"<familyA>_<familyB>_fused" label isn't one, so a fused-ORF hit currently
ships with no calibrated cutoff of its own in <release>.diamond_cutoffs.tsv.
There's no train/test split for fusion candidates to calibrate against
either (they're a small side pool pulled out of 01c, never through
03_split_train_test.py). Whatever osmotool does with an uncalibrated label
is untested. If this turns out to matter in practice, the two obvious
fixes are (a) extend 08b to also process each declared fusion pair's
merged label, using the fusion_candidates pool directly instead of a
train/test split, or (b) have osmotool fall back to the more permissive
(lower) of the pair's two component cutoffs for a "*_fused" hit -- neither
implemented here.

Usage:
  python 08d_build_fusion_refs.py --refs refs --families families.yaml
Output:
  refs/<familyA>_<familyB>_fused.faa   (one per declared, non-empty pair)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def load_fusion_pairs(path: Path) -> list[tuple[str, str]]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    fam_by_name = {fam["name"]: fam for fam in data["families"]}
    seen = set()
    pairs = []
    for fam in data["families"]:
        partner = fam.get("fusion_partner")
        if not partner or partner not in fam_by_name:
            continue
        key = frozenset((fam["name"], partner))
        if key in seen:
            continue
        seen.add(key)
        pairs.append(tuple(sorted((fam["name"], partner))))
    return pairs


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


def header_accession(header: str) -> str:
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header


def retag_as_fused(header: str, label: str) -> str:
    parts = header.split("|", 1)
    rest = parts[1] if len(parts) > 1 else ""
    return f"{label}|{rest}" if rest else label


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    args = ap.parse_args()

    pairs = load_fusion_pairs(args.families)
    if not pairs:
        print("No families.yaml fusion_partner pairs declared -- nothing to build.")
        return

    for fam_a, fam_b in pairs:
        label = f"{fam_a}_{fam_b}_fused"
        by_accession: dict[str, tuple[str, str]] = {}
        for fam in (fam_a, fam_b):
            candidates_path = args.refs / f"{fam}.positive.fusion_candidates.faa"
            for header, seq in parse_fasta(candidates_path):
                by_accession.setdefault(header_accession(header), (header, seq))

        if not by_accession:
            print(f"[{fam_a}/{fam_b}] no fusion candidates found -- skipping "
                  f"(expected if this lineage/pair has no fused-ORF representatives "
                  f"in the fetched data)")
            continue

        retagged = [(retag_as_fused(header, label), seq) for header, seq in by_accession.values()]
        out_path = args.refs / f"{label}.faa"
        write_fasta(out_path, retagged)
        print(f"[{fam_a}/{fam_b}] {len(retagged)} fused-ORF references -> {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
