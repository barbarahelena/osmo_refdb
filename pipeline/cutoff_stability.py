#!/usr/bin/env python3
"""
cutoff_stability.py -- how much does a family's GA cutoff move when the
calibration data change a little, for each way of choosing it?

A rebuild never sees exactly the same calibration sequences (a different
split, a changed positive set), so a cutoff that jumps between nearly
identical calibration sets is not really determined by the data. This repeatedly
draws 80% of a family's calibration positives and negatives without
replacement (hmms/scores/<family>.{positive,negative}.tblout, plus the
decoy contest for decoy-aware families), chooses the cutoff on each draw by
  sweep   the plain F1 sweep (the original method), and
  bagged  the median of bootstrap F1 optima (06_calibrate_cutoffs.py),
and reports, per family, how far each moves across the draws (the
interquartile range, in bits) and what each costs in F1 on the *full*
calibration data compared with the best the full data allow.

Only families whose cutoff comes from an F1 sweep are compared (the others
are set by clean separation or by the largest negative).

Usage:
  python pipeline/cutoff_stability.py --release releases/v10 [--draws 100] [--out stability.tsv]
  python pipeline/cutoff_stability.py --release releases/v10 --compare releases/v10-taxsplit
    (how far each method's cutoff moves between two real builds)
Deterministic: every random choice is seeded from the family name.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import random
import statistics as st
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("calibrate", HERE / "06_calibrate_cutoffs.py")
cal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cal)

ARGS: dict = {}


def work(family: str) -> dict | None:
    score_dir = ARGS["release"] / "hmms" / "scores"
    aware = ARGS["aware"].get(family)
    n_lost = 0
    if aware:
        pos, neg, n_lost, _ = cal.apply_decoys(score_dir, family, aware)
    else:
        pos = cal.parse_tblout_scores(score_dir / f"{family}.positive.tblout")
        neg = cal.parse_tblout_scores(score_dir / f"{family}.negative.tblout")
    if len(pos) == 0 or len(neg) < cal.MIN_NEG_FOR_RELIABLE_F1 or neg.max() < pos.min():
        return None                                   # not an F1-swept family
    full_sweep, full_f1 = cal.choose_cutoff_by_f1(pos, neg, n_lost)
    rng = random.Random(zlib.crc32(f"{family}|stability".encode()))
    d = ARGS["draws"]
    out = {"sweep": [], "bagged": []}
    for _ in range(d):
        p = np.array(rng.sample(list(pos), max(2, int(0.8 * len(pos)))))
        n = np.array(rng.sample(list(neg), max(2, int(0.8 * len(neg)))))
        lost = round(0.8 * n_lost)
        out["sweep"].append(cal.choose_cutoff_by_f1(p, n, lost)[0])
        out["bagged"].append(cal.bagged_cutoff(p, n, lost, family=f"{family}|{_}", n_boot=ARGS["boot"])[0])
    row = {"family": family, "n_pos": len(pos), "n_neg": len(neg), "full_data_best_f1": round(full_f1, 3)}
    for m, cuts in out.items():
        q1, med, q3 = np.percentile(cuts, [25, 50, 75])
        f1s = [cal.f1_at(pos, neg, c, n_lost) for c in cuts]
        row[f"{m}_median"] = round(float(med), 1)
        row[f"{m}_iqr_bits"] = round(float(q3 - q1), 1)
        row[f"{m}_mean_f1_cost"] = round(full_f1 - float(np.mean(f1s)), 4)
    return row


def release_cutoffs(release: Path, family: str, aware: dict) -> dict | None:
    """Both methods' cutoffs for a family's real calibration scores in a release."""
    score_dir = release / "hmms" / "scores"
    n_lost = 0
    if aware.get(family):
        pos, neg, n_lost, _ = cal.apply_decoys(score_dir, family, aware[family])
    else:
        pos = cal.parse_tblout_scores(score_dir / f"{family}.positive.tblout")
        neg = cal.parse_tblout_scores(score_dir / f"{family}.negative.tblout")
    if len(pos) == 0 or len(neg) < cal.MIN_NEG_FOR_RELIABLE_F1 or neg.max() < pos.min():
        return None
    sweep, f1 = cal.choose_cutoff_by_f1(pos, neg, n_lost)
    bagged, p10, p90 = cal.bagged_cutoff(pos, neg, n_lost, family=family)
    return {"n_pos": len(pos), "n_neg": len(neg), "sweep": sweep, "bagged": bagged,
            "p10": p10, "p90": p90, "f1_sweep": f1, "f1_bagged": cal.f1_at(pos, neg, bagged, n_lost)}


def compare(a: Path, b: Path, families_path: Path) -> None:
    """How far each method's cutoff moves between two real builds of the same
    families, whose calibration sets differ (a different split, a changed
    positive set): the quantity a rebuild actually cares about."""
    aware = cal.load_aware_decoys(families_path)
    rows = []
    for fam in cal.load_family_names(families_path):
        x, y = release_cutoffs(a, fam, aware), release_cutoffs(b, fam, aware)
        if x and y:
            rows.append({"family": fam, **{f"{k}_a": v for k, v in x.items()},
                         **{f"{k}_b": v for k, v in y.items()}})
    print(f"cutoffs of {len(rows)} F1-swept families, {a.name} vs {b.name} (same scores recomputed by both methods)\n")
    print(f"{'family':8} {'n_pos a/b':>9} | sweep a -> b  (move) | bagged a -> b  (move) | bagged p10-p90 in a")
    for r in sorted(rows, key=lambda r: -abs(r["sweep_a"] - r["sweep_b"])):
        print(f"{r['family']:8} {r['n_pos_a']:4d}/{r['n_pos_b']:<4d} | {r['sweep_a']:6.1f} -> {r['sweep_b']:6.1f} ({abs(r['sweep_a']-r['sweep_b']):5.1f})"
              f" | {r['bagged_a']:6.1f} -> {r['bagged_b']:6.1f} ({abs(r['bagged_a']-r['bagged_b']):5.1f}) | {r['p10_a']:.0f}-{r['p90_a']:.0f}")
    ms = [abs(r["sweep_a"] - r["sweep_b"]) for r in rows]
    mb = [abs(r["bagged_a"] - r["bagged_b"]) for r in rows]
    print(f"\nmedian move: sweep {st.median(ms):.1f} bits, bagged {st.median(mb):.1f}; "
          f"mean move: sweep {st.mean(ms):.1f}, bagged {st.mean(mb):.1f}; "
          f"bagged moves less in {sum(y < x for x, y in zip(ms, mb))}, more in {sum(y > x for x, y in zip(ms, mb))} of {len(rows)}")
    print(f"mean F1 on the own calibration data: sweep {st.mean(r['f1_sweep_a'] for r in rows):.3f}, "
          f"bagged {st.mean(r['f1_bagged_a'] for r in rows):.3f} (release {a.name})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compare", type=Path, default=None,
                    help="instead of subsampling, compare this release's cutoffs with --release's")
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--draws", type=int, default=100, help="80%% subsamples per family")
    ap.add_argument("--boot", type=int, default=200, help="bootstrap resamples inside each draw")
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    if args.compare:
        compare(args.release, args.compare, args.families)
        return
    ARGS.update(release=args.release, draws=args.draws, boot=args.boot,
                aware=cal.load_aware_decoys(args.families))
    families = cal.load_family_names(args.families)
    with Pool(args.processes) as pool:
        rows = [r for r in pool.map(work, families) if r]
    rows.sort(key=lambda r: -r["sweep_iqr_bits"])
    out = args.out or args.release / "hmms" / "cutoff_stability.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} F1-swept families, {args.draws} draws of 80% each\n")
    print(f"{'family':8} {'n_pos':>5} {'n_neg':>5} | cutoff IQR in bits: sweep -> bagged | median cutoff: sweep / bagged | F1 cost: sweep / bagged")
    for r in rows:
        print(f"{r['family']:8} {r['n_pos']:5d} {r['n_neg']:5d} | {r['sweep_iqr_bits']:8.1f} -> {r['bagged_iqr_bits']:6.1f}"
              f"      | {r['sweep_median']:7.1f} / {r['bagged_median']:7.1f} | {r['sweep_mean_f1_cost']:.3f} / {r['bagged_mean_f1_cost']:.3f}")
    better = sum(r["bagged_iqr_bits"] < r["sweep_iqr_bits"] for r in rows)
    print(f"\nbagged is more stable (smaller IQR) in {better} of {len(rows)} families; "
          f"median IQR sweep {st.median(r['sweep_iqr_bits'] for r in rows):.1f} -> bagged "
          f"{st.median(r['bagged_iqr_bits'] for r in rows):.1f} bits; mean F1 cost "
          f"{st.mean(r['sweep_mean_f1_cost'] for r in rows):.4f} (sweep) vs "
          f"{st.mean(r['bagged_mean_f1_cost'] for r in rows):.4f} (bagged)")
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
