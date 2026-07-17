#!/usr/bin/env python3
"""
08b_calibrate_diamond_cutoffs.py — set per-family DIAMOND bit-score cutoffs.

DIAMOND's identity/query-cover/e-value filters are applied uniformly across
every family (a single --min_identity for the whole database), unlike HMM's
per-family calibrated GA cutoff. This is a real limitation, not an
oversight: a flat 80% identity threshold cannot separate a true family
member from a genuinely close paralog whose overall sequence identity to
it happens to exceed 80%. Found in production: real E. coli BetT/CaiT --
deliberately curated as betL's hard negatives, specifically because
they're choline/carnitine-specific rather than glycine-betaine-specific --
scored 95.6% identity to a betL reference and passed straight through
DIAMOND's flat threshold when nothing analogous to HMM's GA cutoff existed
to catch it.

This calibrates a per-family minimum DIAMOND bitscore, mirroring how
06_calibrate_cutoffs.py calibrates HMM's GA cutoffs: score each family's
held-out positive TEST set and (purity+length-filtered) hard-negative set
against the built DIAMOND database with permissive settings (--id 0
--query-cover 0, so the true score distribution isn't clipped by
production-default filters), restricted to each query's best hit
specifically against THAT family's own references (not just any hit in
the combined db), and pick a cutoff that separates the two distributions
-- same clean_separation / overlapping_distributions_review_needed logic
as the HMM calibration, just in DIAMOND bitscore space.

Usage:
  python 08b_calibrate_diamond_cutoffs.py --refs refs --release releases/v3 \\
      --release-name osmo_refdb --families families.yaml
Output:
  <release>/<release-name>.diamond_cutoffs.tsv
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import yaml

NEG_PERCENTILE = 99.0


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def best_family_restricted_scores(
    query_faa: Path, db: Path, target_family: str, threads: int, tmp_dir: Path
) -> np.ndarray:
    """
    Query *query_faa* against the full combined DIAMOND *db* with permissive
    settings, and return, per query sequence, its best bitscore among hits
    to *target_family*'s own references specifically (0.0 if a query has no
    hit to that family at all, e.g. it matched only some other family or
    nothing). --max-target-seqs 0 (unlimited) so a target_family hit is
    never silently excluded by hits from other, stronger-scoring families.
    """
    if not query_faa.exists() or query_faa.stat().st_size == 0:
        return np.array([])

    out_path = tmp_dir / f"{target_family}.diamond_calib.tsv"
    subprocess.run(
        ["diamond", "blastp",
         "--query", str(query_faa), "--db", str(db),
         "--out", str(out_path), "--outfmt", "6", "qseqid", "sseqid", "bitscore",
         "--id", "0", "--query-cover", "0", "--evalue", "10",
         "--max-target-seqs", "0", "--threads", str(threads), "--quiet"],
        check=True,
    )

    best_per_query: dict[str, float] = {}
    n_queries = 0
    with open(query_faa) as fh:
        for line in fh:
            if line.startswith(">"):
                n_queries += 1
                best_per_query[line[1:].strip()] = 0.0

    with open(out_path) as fh:
        for line in fh:
            qseqid, sseqid, bitscore = line.rstrip("\n").split("\t")
            if sseqid.split("|")[0] != target_family:
                continue
            bitscore = float(bitscore)
            if qseqid in best_per_query and bitscore > best_per_query[qseqid]:
                best_per_query[qseqid] = bitscore

    return np.array(list(best_per_query.values()))


def choose_cutoff(pos_scores: np.ndarray, neg_scores: np.ndarray) -> tuple[float, str]:
    """Same logic as 06_calibrate_cutoffs.py's choose_cutoff(), applied to
    DIAMOND bitscores instead of HMM bitscores."""
    if len(pos_scores) == 0:
        return 0.0, "no_positive_scores"

    min_pos = float(np.min(pos_scores))

    if len(neg_scores) == 0:
        return max(min_pos - 1.0, 0.0), "no_negatives_available"

    max_neg = float(np.max(neg_scores))

    if max_neg < min_pos:
        cutoff = (max_neg + min_pos) / 2.0
        return cutoff, "clean_separation"
    else:
        cutoff = float(np.percentile(neg_scores, NEG_PERCENTILE))
        return cutoff, "overlapping_distributions_review_needed"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--release", type=Path, required=True,
                     help="Release directory containing <release-name>.dmnd")
    ap.add_argument("--release-name", default="osmo_refdb")
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    db = args.release / f"{args.release_name}.dmnd"
    if not db.exists():
        print(f"ERROR: {db} not found -- run 08_build_diamond_db.sh first.")
        return

    families = load_family_names(args.families)
    manifest_rows = []

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for family in families:
            pos_test = args.refs / f"{family}.positive.test.faa"
            neg = args.refs / f"{family}.negative.faa"

            pos_scores = best_family_restricted_scores(pos_test, db, family, args.threads, tmp_dir)
            neg_scores = best_family_restricted_scores(neg, db, family, args.threads, tmp_dir)

            cutoff, status = choose_cutoff(pos_scores, neg_scores)

            print(f"[{family}] cutoff={cutoff:.1f} status={status} "
                  f"(n_pos={len(pos_scores)}, n_neg={len(neg_scores)})")

            manifest_rows.append({
                "family": family,
                "cutoff_bitscore": round(cutoff, 1),
                "status": status,
                "n_positive_scored": len(pos_scores),
                "n_negative_scored": len(neg_scores),
                "min_positive_score": round(float(np.min(pos_scores)), 1) if len(pos_scores) else "",
                "max_negative_score": round(float(np.max(neg_scores)), 1) if len(neg_scores) else "",
            })

    if not manifest_rows:
        print("No families calibrated.")
        return

    manifest_path = args.release / f"{args.release_name}.diamond_cutoffs.tsv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    n_review = sum(1 for r in manifest_rows if r["status"] not in ("clean_separation",))
    print(f"\nDone. {len(manifest_rows)} families calibrated, {n_review} flagged for review.")
    print(f"Manifest: {manifest_path}")
    print("Ship this file alongside the .dmnd release; osmotool consumes it via "
          "--diamond_cutoffs.")


if __name__ == "__main__":
    main()
