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
  * if there IS overlap: sweep every candidate threshold across the
    combined positive+negative score population and pick the one that
    maximizes F1 (same technique 11_compute_metrics.py's
    summarize_best_threshold() already uses for post-hoc benchmark
    analysis, applied here at calibration time instead).

    This replaces an earlier flat "99th percentile of negatives" rule,
    which ignored the positive distribution entirely: for a family whose
    surviving hard negatives remain close to true positives even after
    01b/01c's purity+length QC (e.g. mscL, where 72% of raw negatives
    were flagged as contaminated and what survived was still close
    enough to true mscL to push a percentile-of-negatives cutoff high),
    that rule could reject a large fraction of genuine positives just to
    protect against an already-mostly-cleaned negative pool. Confirmed in
    production: mscL's old cutoff rejected a real, 100%-identical protein
    outright when tested against a real genome (Rhodopirellula baltica).
    The F1 sweep balances both distributions instead of only bounding
    false positives, and this family is still flagged in the manifest for
    manual review — the benchmark (ROC/PR) remains the final word on
    whether the family is usable with an HMM.

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


MIN_NEG_FOR_RELIABLE_F1 = 30


def choose_cutoff(pos_scores: np.ndarray, neg_scores: np.ndarray, neg_file_exists: bool) -> tuple[float, str, float | None]:
    if len(pos_scores) == 0:
        return 0.0, "no_positive_scores", None

    min_pos = float(np.min(pos_scores))

    if len(neg_scores) == 0:
        if neg_file_exists:
            # Negatives were fetched and searched, but none scored above
            # HMMER's default reporting threshold (E<=10) — i.e. the HMM is
            # highly specific and doesn't detect the negative family at all.
            # This is the best possible outcome, not a missing-data case.
            return max(min_pos - 1.0, 0.0), "clean_separation_no_negative_hits", 1.0
        return max(min_pos - 1.0, 0.0), "no_negatives_available", None

    max_neg = float(np.max(neg_scores))

    if max_neg < min_pos:
        cutoff = (max_neg + min_pos) / 2.0
        return cutoff, "clean_separation", 1.0
    elif len(neg_scores) < MIN_NEG_FOR_RELIABLE_F1:
        # Too few negatives to trust an F1 sweep. With e.g. only 10 points
        # against hundreds of positives, "F1-optimal" can degenerate to
        # "accept everything" -- a handful of false positives barely dents
        # F1 when true positives vastly outnumber them in the calibration
        # set, but that class balance is an artifact of how few negatives
        # survived curation, not a reflection of the real ratio of
        # "family member" to "everything else" in an actual genome.
        # Confirmed in production: proP's 10-negative F1 sweep produced a
        # cutoff that accepted every calibration positive, then let
        # through several weak, high-compositional-bias false positives
        # on a real genome. Fall back to requiring the cutoff clear every
        # observed negative instead -- conservative, and flagged as a
        # data problem (see families.yaml's negative_query design notes
        # for this family) rather than papered over with a cleverer
        # formula, since no cutoff-selection algorithm substitutes for a
        # negative set too small to be representative.
        cutoff = max_neg + 1.0
        return cutoff, "insufficient_negative_data", None
    else:
        cutoff, f1 = choose_cutoff_by_f1(pos_scores, neg_scores)
        return cutoff, "overlapping_distributions_f1_calibrated", f1


def choose_cutoff_by_f1(pos_scores: np.ndarray, neg_scores: np.ndarray) -> tuple[float, float]:
    """
    Sweep every distinct score in the combined positive+negative
    population as a candidate cutoff (>=) and return the one that
    maximizes F1, plus the F1 achieved at that cutoff.

    Ties (multiple thresholds achieving the same best F1) are broken by
    preferring the LOWEST such threshold, i.e. the most permissive one
    that still achieves the best balance -- favors recall when the
    balance point is genuinely ambiguous, rather than an arbitrary
    numpy-sort-order tiebreak.
    """
    n_pos_total = len(pos_scores)
    candidates = sorted(set(pos_scores.tolist()) | set(neg_scores.tolist()))

    best_f1 = -1.0
    best_thr = candidates[0] if candidates else 0.0
    for thr in candidates:
        tp = int((pos_scores >= thr).sum())
        fn = n_pos_total - tp
        fp = int((neg_scores >= thr).sum())
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
        if f1 > best_f1:
            best_f1 = f1
            best_thr = thr

    return best_thr, best_f1


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
        cutoff, status, f1_at_cutoff = choose_cutoff(pos_scores, neg_scores, neg_file_exists)

        set_ga_line(hmm_path, cutoff)

        f1_str = f"{f1_at_cutoff:.3f}" if f1_at_cutoff is not None else "n/a"
        print(f"[{family}] cutoff={cutoff:.2f} status={status} f1={f1_str} "
              f"(n_pos={len(pos_scores)}, n_neg={len(neg_scores)})")

        manifest_rows.append({
            "family": family,
            "cutoff_bits": round(cutoff, 2),
            "status": status,
            "f1_at_cutoff": round(f1_at_cutoff, 3) if f1_at_cutoff is not None else "",
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
