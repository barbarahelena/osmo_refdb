#!/usr/bin/env python3
"""
06b_qc_scorecard.py — consolidate per-family QC signal into one scorecard.

Three separate checks each write their own manifest:
  - 01b_check_negative_purity.py  -> refs/negative_purity_manifest.tsv
  - 01c_check_length_outliers.py  -> refs/length_outlier_manifest.tsv
  - 06_calibrate_cutoffs.py       -> hmms/cutoff_manifest.tsv
Cross-referencing them by hand (as done manually earlier in this project,
e.g. to diagnose proX/betL) works but doesn't scale as more families get
added. This script merges all three into one per-family table so you can
see at a glance which families are solid vs. still need curation, instead
of opening three files.

A family is flagged "review_needed" if any of:
  - >10% of its hard negatives were flagged as likely mislabeled orthologs
  - >10% of its positive OR negative sequences were flagged as length
    outliers (likely fusion proteins/fragments)
  - its HMM cutoff status isn't a clean-separation status

Usage:
  python 06b_qc_scorecard.py --refs refs --hmms hmms --families families.yaml \\
      --out qc_scorecard.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml

FLAG_FRACTION_THRESHOLD = 0.10
CLEAN_HMM_STATUSES = {"clean_separation", "clean_separation_no_negative_hits"}


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def load_tsv_by_family(path: Path, key_field: str = "family") -> dict[str, list[dict]]:
    """Load a manifest TSV into {family: [row, ...]} (list, since
    length_outlier_manifest.tsv has one row per family+label)."""
    rows_by_family: dict[str, list[dict]] = {}
    if not path.exists():
        return rows_by_family
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows_by_family.setdefault(row[key_field], []).append(row)
    return rows_by_family


def frac(numerator: str, denominator: str) -> float:
    n, d = int(numerator), int(denominator)
    return (n / d) if d else 0.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--hmms", type=Path, default=Path("hmms"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--out", type=Path, default=Path("qc_scorecard.tsv"))
    args = ap.parse_args()

    families = load_family_names(args.families)
    purity = load_tsv_by_family(args.refs / "negative_purity_manifest.tsv")
    length = load_tsv_by_family(args.refs / "length_outlier_manifest.tsv")
    cutoffs = load_tsv_by_family(args.hmms / "cutoff_manifest.tsv")

    if not purity and not length and not cutoffs:
        print("No source manifests found (negative_purity_manifest.tsv, "
              "length_outlier_manifest.tsv, cutoff_manifest.tsv) -- run "
              "01b/01c/06 first.")
        return

    rows = []
    for family in families:
        p = purity.get(family, [{}])[0]
        length_rows = {r.get("label"): r for r in length.get(family, [])}
        pos_len = length_rows.get("positive", {})
        neg_len = length_rows.get("negative", {})
        c = cutoffs.get(family, [{}])[0]

        neg_purity_frac = frac(p["n_flagged"], p["n_total_negatives"]) if p else None
        pos_length_frac = frac(pos_len["n_flagged"], pos_len["n_total"]) if pos_len else None
        neg_length_frac = frac(neg_len["n_flagged"], neg_len["n_total"]) if neg_len else None
        hmm_status = c.get("status", "")

        review_reasons = []
        if neg_purity_frac is not None and neg_purity_frac > FLAG_FRACTION_THRESHOLD:
            review_reasons.append(f"neg_purity={neg_purity_frac:.0%}")
        if pos_length_frac is not None and pos_length_frac > FLAG_FRACTION_THRESHOLD:
            review_reasons.append(f"pos_length={pos_length_frac:.0%}")
        if neg_length_frac is not None and neg_length_frac > FLAG_FRACTION_THRESHOLD:
            review_reasons.append(f"neg_length={neg_length_frac:.0%}")
        if hmm_status and hmm_status not in CLEAN_HMM_STATUSES:
            review_reasons.append(f"hmm_status={hmm_status}")

        rows.append({
            "family": family,
            "neg_purity_flagged": p.get("n_flagged", ""),
            "neg_purity_total": p.get("n_total_negatives", ""),
            "pos_length_flagged": pos_len.get("n_flagged", ""),
            "pos_length_total": pos_len.get("n_total", ""),
            "neg_length_flagged": neg_len.get("n_flagged", ""),
            "neg_length_total": neg_len.get("n_total", ""),
            "hmm_cutoff_bits": c.get("cutoff_bits", ""),
            "hmm_cutoff_status": hmm_status,
            "review_needed": "yes" if review_reasons else "no",
            "review_reasons": ";".join(review_reasons),
        })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_review = sum(1 for r in rows if r["review_needed"] == "yes")
    print(f"{len(rows)} families scored, {n_review} flagged for review.")
    for r in rows:
        if r["review_needed"] == "yes":
            print(f"  [{r['family']}] {r['review_reasons']}")
    print(f"\nScorecard: {args.out}")


if __name__ == "__main__":
    main()
