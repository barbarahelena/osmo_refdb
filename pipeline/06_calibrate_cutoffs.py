#!/usr/bin/env python3
"""
06_calibrate_cutoffs.py — set per-family HMM gathering thresholds (GA cutoffs)

Reads the hmmsearch --tblout score tables produced by 05_build_hmms.sh for
the positive calibration set (positive.calib: held out from the model and
from the benchmark's positive.test) and hard-negative set, finds a score
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

families.yaml `pfam_model: PFxxxxx` families (see 05_build_hmms.sh) are
handled differently here: their .hmm file already carries a real GA cutoff
-- Pfam's own, curated for that Pfam family, not computed by this script --
so it's read and left as-is rather than overwritten. Our own held-out
positive/negative sets are still scored against it (already done by
05_build_hmms.sh) purely to confirm Pfam's cutoff actually separates OUR
curated true positives from OUR hard negatives; a bad result there
(status="pfam_ga_review_needed" in the manifest) means the assumption that
this Pfam family is gene-specific enough for our purposes was wrong, not
that the number itself needs recalculating -- if that happens, remove
pfam_model from that family in families.yaml and let it fall back to a
normal custom-built, empirically-calibrated model instead.

Decoy-aware GA: for a family a `decoy_aware_ga` HMM decoy competes with
(families.yaml hmm_decoys; scored by 05b_build_hmm_decoys.py), a calibration
sequence that the decoy outscores (and that reaches the decoy's own GA) is
not called as the family, as in annotate's best-hit contest: such a
positive counts as a miss at every cutoff, such a negative drops out. The
cutoff then only has to separate what the decoy leaves, which lets a family
whose confusers overlap it in score (betL vs betT/caiT) keep a lower GA.

Usage: python 06_calibrate_cutoffs.py [--hmms hmms] [--families families.yaml]
Output: updates hmms/<family>.hmm in place with a GA line (skipped for
        pfam_model families, whose GA line is Pfam's own)
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


def load_pfam_model_families(path: Path) -> dict[str, str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return {fam["name"]: fam["pfam_model"] for fam in data["families"] if fam.get("pfam_model")}


def get_ga_line(hmm_path: Path) -> float | None:
    match = re.search(r"^GA\s+([\d.]+)", hmm_path.read_text(), flags=re.MULTILINE)
    return float(match.group(1)) if match else None


def parse_tblout_named(path: Path) -> dict[str, float]:
    """{sequence name: full-sequence bit score} from an hmmsearch --tblout."""
    out: dict[str, float] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            out[f[0]] = max(out.get(f[0], float("-inf")), float(f[5]))
    return out


def load_aware_decoys(path: Path) -> dict[str, list[tuple[str, float]]]:
    """{family: [(decoy name, decoy GA)]} for hmm_decoys with decoy_aware_ga."""
    with open(path) as fh:
        data = yaml.safe_load(fh)
    out: dict[str, list[tuple[str, float]]] = {}
    for d in data.get("hmm_decoys") or []:
        if d.get("decoy_aware_ga"):
            for fam in d.get("competes_with", []):
                out.setdefault(fam, []).append((d["name"], float(d["ga"])))
    return out


def apply_decoys(score_dir: Path, family: str, decoys: list[tuple[str, float]]
                 ) -> tuple[np.ndarray, np.ndarray, int, int]:
    """The family's calibration scores after the decoy contest: (positive
    scores the family keeps, negative scores it keeps, positives lost to a
    decoy, negatives won by a decoy)."""
    kept = {}
    lost = {}
    for kind in ("positive", "negative"):
        own = parse_tblout_named(score_dir / f"{family}.{kind}.tblout")
        rival: dict[str, float] = {}
        for name, ga in decoys:
            for seq, sc in parse_tblout_named(score_dir / f"{family}.{name}.{kind}.tblout").items():
                if sc >= ga:
                    rival[seq] = max(rival.get(seq, float("-inf")), sc)
        won = {s for s, sc in own.items() if rival.get(s, float("-inf")) > sc}
        kept[kind] = np.array([sc for s, sc in own.items() if s not in won])
        lost[kind] = len(won)
    return kept["positive"], kept["negative"], lost["positive"], lost["negative"]


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


def choose_cutoff(pos_scores: np.ndarray, neg_scores: np.ndarray, neg_file_exists: bool,
                  n_pos_lost: int = 0) -> tuple[float, str, float | None]:
    """n_pos_lost: positives a decoy takes (decoy-aware GA); they count as
    misses at every cutoff, so they lower F1 but never move the cutoff."""
    if n_pos_lost:
        if len(pos_scores) == 0:
            return 0.0, "no_positive_scores", None
        if len(neg_scores) == 0:
            cutoff = max(float(np.min(pos_scores)) - 1.0, 0.0)
            return cutoff, "decoy_aware_clean", 2 * len(pos_scores) / (2 * len(pos_scores) + n_pos_lost)
        cutoff, f1 = choose_cutoff_by_f1(pos_scores, neg_scores, n_pos_lost)
        return cutoff, "decoy_aware_f1_calibrated", f1
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


def choose_cutoff_by_f1(pos_scores: np.ndarray, neg_scores: np.ndarray,
                        n_pos_lost: int = 0) -> tuple[float, float]:
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
    n_pos_total = len(pos_scores) + n_pos_lost
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
    pfam_model_families = load_pfam_model_families(args.families)
    aware_decoys = load_aware_decoys(args.families)

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
        n_pos_lost = n_neg_won = 0
        if family in aware_decoys:
            pos_scores, neg_scores, n_pos_lost, n_neg_won = apply_decoys(
                score_dir, family, aware_decoys[family])
            print(f"[{family}] decoy-aware: {', '.join(n for n, _ in aware_decoys[family])} "
                  f"takes {n_pos_lost} calibration positives and {n_neg_won} negatives")

        if family in pfam_model_families:
            cutoff = get_ga_line(hmm_path)
            if cutoff is None:
                print(f"[{family}] WARNING: pfam_model={pfam_model_families[family]} but no GA line "
                      f"found in {hmm_path} -- falling back to empirical calibration")
                cutoff, status, f1_at_cutoff = choose_cutoff(pos_scores, neg_scores, neg_file_exists)
                set_ga_line(hmm_path, cutoff)
            else:
                # Don't overwrite -- this IS Pfam's own GA line already.
                # Just check how well it separates OUR curated sets, for
                # the manifest (see module docstring).
                n_pos_below = int((pos_scores < cutoff).sum()) if len(pos_scores) else 0
                n_neg_above = int((neg_scores >= cutoff).sum()) if len(neg_scores) else 0
                status = ("pfam_ga_clean" if n_pos_below == 0 and n_neg_above == 0
                          else "pfam_ga_review_needed")
                tp = len(pos_scores) - n_pos_below
                fp = n_neg_above
                fn = n_pos_below
                precision = tp / (tp + fp) if (tp + fp) else 0.0
                recall = tp / (tp + fn) if (tp + fn) else 0.0
                f1_at_cutoff = (2 * precision * recall / (precision + recall)
                                 if (precision + recall) else None)
        else:
            cutoff, status, f1_at_cutoff = choose_cutoff(pos_scores, neg_scores, neg_file_exists,
                                                         n_pos_lost)
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
            "positives_lost_to_decoy": n_pos_lost,
            "negatives_won_by_decoy": n_neg_won,
        })

    manifest_path = hmm_dir / "cutoff_manifest.tsv"
    if not manifest_rows:
        print("No families calibrated (no .hmm files found under hmms/).")
        return
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    CLEAN_STATUSES = {"clean_separation", "pfam_ga_clean"}
    n_review = sum(1 for r in manifest_rows if r["status"] not in CLEAN_STATUSES)
    print(f"\nDone. {len(manifest_rows)} families calibrated, {n_review} flagged for review.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
