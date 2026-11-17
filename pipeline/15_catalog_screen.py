#!/usr/bin/env python3
"""
15_catalog_screen.py -- screen a release's HMMs for confuser classes by
annotating a real gene catalog, before the release is frozen.

Why: the benchmark's hard negatives only contain the confusers someone
thought of (gor's are trxB only), so a family can look fine there and still
claim thousands of genes that are not its own. A large real-world catalog is
the widest specificity test there is, and it is run anyway.

For every panel family it reports how many catalog genes it claimed and, the
useful part, how many of those score *below anything one of the family's own
true positives scored* in calibration (hmms/scores/<family>.positive.tblout,
the held-out positive.calib set). A family whose claimed genes mostly score
below its weakest true positive is calling something else: the cutoff sits
in a gap that true members never occupy. Example, v10 on 159,560 real human
gut genes: gor claimed 127 genes (the largest family, 18% of all 691
annotated genes), 86% of them below its weakest true positive (404 bits,
cutoff 245.6); every other family with >= 20 genes had <= 5% below.

Output: <catalog>/screen.tsv, a printed table, and the number of genes each
HMM decoy took (a decoy that takes nothing is not needed; one that takes a lot
is doing its job). A flag is a prompt to look, not a verdict: families with
few genes are listed but not flagged, and a broad-by-design family (cspA) can
score high on the looser measure.

Usage:
  python pipeline/15_catalog_screen.py --release releases/v10 --catalog gmgc_gut
  python pipeline/15_catalog_screen.py --release releases/v10 --catalog-dir /path/to/catalog_dir
"""
from __future__ import annotations

import argparse
import csv
import statistics as st
from pathlib import Path

import yaml


def calib_scores(path: Path) -> list[float]:
    if not path.exists():
        return []
    return sorted(float(l.split()[5]) for l in path.read_text().splitlines()
                  if l.strip() and not l.startswith("#"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--catalog", help="name of a catalog under <release>/catalogs/")
    ap.add_argument("--catalog-dir", type=Path, help="a catalog directory (instead of --catalog)")
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--min-genes", type=int, default=20,
                    help="families with fewer claimed genes are listed but never flagged")
    ap.add_argument("--flag-share", type=float, default=0.25,
                    help="flag a family when this share of its genes scores below its "
                         "weakest calibration positive")
    args = ap.parse_args()
    cdir = args.catalog_dir or args.release / "catalogs" / (args.catalog or "")
    annotation = cdir / "catalog.annotation.tsv"
    if not annotation.exists():
        raise SystemExit(f"ERROR: {annotation} not found (use --catalog NAME or --catalog-dir)")

    decoys = {d["name"]: d.get("competes_with", [])
              for d in (yaml.safe_load(open(args.families)).get("hmm_decoys") or [])} \
        if args.families.exists() else {}
    claimed: dict[str, list[float]] = {}
    taken: dict[str, int] = {}
    n_genes = n_hit = n_homolog = 0
    with open(annotation) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_hit += 1
            if r["status"] == "annotated":
                claimed.setdefault(r["label"], []).append(float(r["score"]))
            elif r["status"] == "decoy":
                taken[r["label"]] = taken.get(r["label"], 0) + 1
            else:
                n_homolog += 1
    n_annot = sum(len(v) for v in claimed.values())

    rows = []
    for fam, sc in claimed.items():
        pos = calib_scores(args.release / "hmms" / "scores" / f"{fam}.positive.tblout")
        if len(pos) < 10:
            rows.append({"family": fam, "genes": len(sc), "share_of_annotated": round(len(sc) / n_annot, 4),
                         "pct_below_weakest_true_positive": "", "pct_below_true_positive_p5": "",
                         "median_score_claimed": round(st.median(sc), 1), "weakest_true_positive": "",
                         "flag": "no calibration positives"})
            continue
        p5 = pos[int(0.05 * len(pos))]
        below_min = sum(s < pos[0] for s in sc) / len(sc)
        below_p5 = sum(s < p5 for s in sc) / len(sc)
        flag = ("CHECK: likely confuser class" if len(sc) >= args.min_genes and below_min >= args.flag_share
                else "too few genes" if len(sc) < args.min_genes else "")
        rows.append({"family": fam, "genes": len(sc), "share_of_annotated": round(len(sc) / n_annot, 4),
                     "pct_below_weakest_true_positive": round(100 * below_min, 1),
                     "pct_below_true_positive_p5": round(100 * below_p5, 1),
                     "median_score_claimed": round(st.median(sc), 1),
                     "weakest_true_positive": round(pos[0], 1), "flag": flag})
    rows.sort(key=lambda r: -r["genes"])
    out = cdir / "screen.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{annotation.parent.name}: {n_hit} genes with an HMM hit: {n_annot} claimed by a panel family, "
          f"{n_homolog} homologs below cutoff, {sum(taken.values())} taken by decoys")
    print(f"\n{'family':8} {'genes':>7} {'share':>6} | % below weakest true positive / below p5 | median claimed | flag")
    for r in rows[:25]:
        b = f"{r['pct_below_weakest_true_positive']:>5}% / {r['pct_below_true_positive_p5']:>5}%" \
            if r["pct_below_weakest_true_positive"] != "" else "      n/a"
        print(f"{r['family']:8} {r['genes']:7d} {100*r['share_of_annotated']:5.1f}% | {b:>30} | "
              f"{r['median_score_claimed']:>7} | {r['flag']}")
    if taken or decoys:
        print("\nHMM decoys: genes taken (a decoy that takes nothing is not needed)")
        for name in sorted(set(decoys) | set(taken)):
            print(f"  {name:22} {taken.get(name, 0):7d}   competes with {', '.join(decoys.get(name, [])) or '-'}")
    flagged = [r["family"] for r in rows if r["flag"].startswith("CHECK")]
    print(f"\nFlagged: {', '.join(flagged) if flagged else 'none'}   -> {out}")
    print("A flag is a prompt to look (what are those genes?), not a verdict.")


if __name__ == "__main__":
    main()
