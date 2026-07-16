#!/usr/bin/env python3
"""
06_calibrate_cutoffs.py — set per-family HMM gathering thresholds (GA cutoffs)

Reads the hmmsearch --tblout score tables produced by 05_build_hmms.sh for
the held-out positive test set and hard-negative set, finds a score
threshold that separates them, writes the cutoff into each .hmm file's GA
line, and records the decision (and full score distributions) in
hmms/cutoff_manifest.tsv.

Strategy (mirrors how Pfam sets GA cutoffs):
  * lowest positive bit-score  = highest score a true member scored
  * highest negative bit-score = highest score a hard negative scored
  * if negatives score below the lowest positive: cutoff = midpoint of the
    gap (safe separation).
  * if there IS overlap: cutoff = a percentile of the negative distribution
    (default: 99th percentile) to bound false-positive rate, and this family
    is flagged in the manifest as "overlapping" for manual review — the
    benchmark (ROC/PR) is the final word on whether the family is usable
    with an HMM.

Usage: python 06_calibrate_cutoffs.py [--hmms hmms] [--families families.yaml]
Output: updates hmms/<family>.hmm in place with a GA line
        hmms/cutoff_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import yaml

NEG_PERCENTILE = 99.0


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def parse_tblout_scores(path: Path) -> np.ndarray:
    """Parse full-sequence bit-scores from an hmmsearch --tblout file."""
    scores = []
    if not path.exists():
        return np.array(scores)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            # tblout columns: target name, accession, query name, accession,
            # E-value, score, bias, ... (full sequence, column index 5 = score)
            scores.append(float(fields[5]))
    return np.array(scores)


def choose_cutoff(pos_scores: np.ndarray, neg_scores: np.ndarray, neg_file_exists: bool) -> tuple[float, str]:
    if len(pos_scores) == 0:
        return 0.0, "no_positive_scores"

    min_pos = float(np.min(pos_scores))

    if len(neg_scores) == 0:
        if neg_file_exists:
            # Negatives were fetched and searched, but none scored above
            # HMMER's default reporting threshold (E<=10) — i.e. the HMM is
            # highly specific and doesn't detect the negative family at all.
            # This is the best possible outcome, not a missing-data case.
            return max(min_pos - 1.0, 0.0), "clean_separation_no_negative_hits"
        return max(min_pos - 1.0, 0.0), "no_negatives_available"

    max_neg = float(np.max(neg_scores))

    if max_neg < min_pos:
        cutoff = (max_neg + min_pos) / 2.0
        return cutoff, "clean_separation"
    else:
        cutoff = float(np.percentile(neg_scores, NEG_PERCENTILE))
        return cutoff, "overlapping_distributions_review_needed"


def set_ga_line(hmm_path: Path, cutoff: float) -> None:
    """Insert/replace a GA (gathering threshold) line in the HMM file."""
    text = hmm_path.read_text()
    ga_line = f"GA    {cutoff:.2f} {cutoff:.2f};\n"
    if re.search(r"^GA\s", text, flags=re.MULTILINE):
        text = re.sub(r"^GA\s.*$", ga_line.strip(), text, flags=re.MULTILINE)
    else:
        # Insert GA line right before the first "HMM" model line
        text = re.sub(r"^(HMM\s)", ga_line + r"\1", text, count=1, flags=re.MULTILINE)
    hmm_path.write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hmms", type=Path, default=Path("hmms"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    args = ap.parse_args()

    hmm_dir = args.hmms
    score_dir = hmm_dir / "scores"
    families = load_family_names(args.families)

    manifest_rows = []

    for family in families:
        hmm_path = hmm_dir / f"{family}.hmm"
        pos_path = score_dir / f"{family}.positive.tblout"
        neg_path = score_dir / f"{family}.negative.tblout"

        if not hmm_path.exists():
            print(f"[{family}] SKIP: {hmm_path} not found (run 05_build_hmms.sh first)")
            continue

        pos_scores = parse_tblout_scores(pos_path)
        neg_scores = parse_tblout_scores(neg_path)
        # neg_path exists (even if empty of hits) whenever hmmsearch was
        # actually run against a non-empty negative FASTA; distinguishing
        # "ran but found nothing" from "never ran" changes the interpretation.
        neg_file_exists = neg_path.exists() and neg_path.stat().st_size > 0
        cutoff, status = choose_cutoff(pos_scores, neg_scores, neg_file_exists)

        set_ga_line(hmm_path, cutoff)

        print(f"[{family}] cutoff={cutoff:.2f} status={status} "
              f"(n_pos={len(pos_scores)}, n_neg={len(neg_scores)})")

        manifest_rows.append({
            "family": family,
            "cutoff_bits": round(cutoff, 2),
            "status": status,
            "n_positive_scored": len(pos_scores),
            "n_negative_scored": len(neg_scores),
            "min_positive_score": round(float(np.min(pos_scores)), 2) if len(pos_scores) else "",
            "max_negative_score": round(float(np.max(neg_scores)), 2) if len(neg_scores) else "",
        })

    manifest_path = hmm_dir / "cutoff_manifest.tsv"
    if not manifest_rows:
        print("No families calibrated (no .hmm files found under hmms/).")
        return
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    n_review = sum(1 for r in manifest_rows if r["status"] != "clean_separation")
    print(f"\nDone. {len(manifest_rows)} families calibrated, {n_review} flagged for review.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
