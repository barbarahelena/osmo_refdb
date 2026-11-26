#!/usr/bin/env python3
"""
compare_indexes.py -- the study's results with two catalog indexes side by side:
the release's shipped (neighbours) index and the full index of every catalog gene,
scored on the same reads (score.py's scores_by_group.tsv in work/<tag>/ for each).

One table per metric (recall, false calls, mates the catalog settled);
rows are the strain groups, columns the method and index.

Usage:
  python compare_indexes.py --work work --small gmgc_gut --full gmgc_full
Writes work/index_comparison.md and prints it.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

GROUPS = [
    ("in", "Strains in the catalogue"),
    ("partial", "Strains partly in the catalogue"),
    ("out_related", "Not in the catalogue, same genus as a catalogued strain"),
    ("out_distant", "Not in the catalogue, distant"),
    ("all", "Whole sample"),
]


def load(path: Path) -> dict[tuple[str, str], dict]:
    with open(path) as fh:
        return {(r["method"], r["group"]): r for r in csv.DictReader(fh, delimiter="\t")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, default=Path(__file__).resolve().parent / "work")
    ap.add_argument("--small", default="gmgc_gut", help="tag of the shipped index's run")
    ap.add_argument("--full", default="gmgc_full", help="tag of the full index's run")
    args = ap.parse_args()

    small = load(args.work / args.small / "scores_by_group.tsv")
    full = load(args.work / args.full / "scores_by_group.tsv")
    cols = [("DIAMOND only (no catalogue)", small, "diamond_only")]
    for method, name in [("current", "current rule"), ("ambiguous_to_diamond", "ambiguous reads to DIAMOND")]:
        cols += [(f"Catalogue, {name}, shipped index", small, method),
                 (f"Catalogue, {name}, full index", full, method)]

    # the DIAMOND-only column must not depend on the index: same reads, same DIAMOND run
    for g, _ in GROUPS:
        a, b = small[("diamond_only", g)], full[("diamond_only", g)]
        if (a["tp"], a["fp"]) != (b["tp"], b["fp"]):
            raise SystemExit(f"DIAMOND only differs between the two runs in group {g}: not the same reads")

    def table(title, cell, use=cols):
        out = [f"**{title}**", "", "| Strains | " + " | ".join(c[0] for c in use) + " |",
               "|---|" + "---:|" * len(use)]
        for g, label in GROUPS:
            out.append(f"| {label} | " + " | ".join(cell(t[(m, g)]) for _, t, m in use) + " |")
        return out + [""]

    lines = [f"# Shipped index ({args.small}) vs full index ({args.full})", ""]
    lines += table("Reads from panel genes that get the right family (recall)",
                   lambda r: f"{float(r['recall']):.2f} ({r['tp']}/{r['panel_mates']})")
    lines += table("False calls (reads given a family they are not from)", lambda r: f"{int(r['fp']):,}")
    # the same mates are settled under both rules (an ambiguous mate is settled, then sent on or not)
    lines += table("Mates the catalogue settled",
                   lambda r: f"{int(r['settled_by_catalog']):,}",
                   use=[("Shipped index", small, "current"), ("Full index", full, "current")])
    text = "\n".join(lines)
    (args.work / "index_comparison.md").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
