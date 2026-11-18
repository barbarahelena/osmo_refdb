#!/usr/bin/env python3
"""
score.py -- false panel calls per million reads, and recall, of `osmotool
profile` on reads from real catalog genes, with the catalog indexed in full
(all), with only its HMM-hit genes (hits), with those plus the genes close to
them (neighbours, find_neighbours.py), and without a catalog (none).

Truth: the source gene's full-length HMM decision from build-catalog
(catalog.annotation.tsv of the slice): a gene "annotated" as a family is
that family (on a fused gene, any family that passes on it); every other
gene (a homolog below the cutoff, a decoy, or no hit at all) is no panel
family. The same annotation decides the catalog's calls, so a catalog call
on a read from its own gene is correct by construction; what the test
measures is where reads go when their gene is (all) or is not (hits) in the
index, and what the translated fallback makes of the rest.

Each read's call is rebuilt per mate: the catalog's decision when the
catalog settled it (call -> its family; drop or ambiguous -> no call), else
its DIAMOND best hit after the cascade (compute_metrics.py's rebuild), no
call for an excluded family (decoys). The rebuild is checked against the
run's gene_counts.tsv.

Usage:
  python score.py --compute-metrics compute_metrics.py --release release/v9 \\
      --annotation catalogs/all/catalog.annotation.tsv --runs runs --reads reads --out results
  (runs/<divergence>/<arm>/, reads/<divergence>.truth.tsv)
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
from collections import Counter, defaultdict
from pathlib import Path

ARMS = ("none", "hits", "neighbours", "all")


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("compute_metrics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def gene_truth(annotation: Path) -> dict[str, frozenset[str]]:
    """{gene: families it counts as}; genes absent are no panel family."""
    truth = {}
    with open(annotation) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["status"] != "annotated":
                continue
            fams = {row["label"]}
            for part in filter(None, row.get("other_envelopes", "").split(",")):
                fams.add(part.split(":")[0])
            truth[row["gene"]] = frozenset(fams)
    return truth


def read_truth(path: Path, genes: dict[str, frozenset[str]]) -> dict[str, frozenset[str]]:
    """{read id with mate (r<k>/1, r<k>/2): families}; empty = no panel family."""
    out = {}
    empty = frozenset()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fams = genes.get(row["gene"], empty)
            out[f"{row['read']}/1"] = fams
            out[f"{row['read']}/2"] = fams
    return out


def load_excluded(release: Path) -> set[str]:
    f = release / "osmo_refdb.profile_excluded_families.txt"
    if not f.exists():
        return set()
    return {line.strip() for line in f.read_text().splitlines() if line.strip() and not line.startswith("#")}


def catalog_calls(run: Path) -> tuple[dict[str, str | None], Counter]:
    """{read/mate: family or None} for every read the catalog settled."""
    settled: dict[str, str | None] = {}
    decisions: Counter = Counter()
    for f in run.glob("*.catalog_hits.tsv"):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rid = f"{row['read']}/{row['mate']}" if row["mate"] in ("1", "2") else row["read"]
                settled[rid] = row["label"] if row["decision"] == "call" else None
                decisions[row["decision"]] += 1
    return settled, decisions


def translated_calls(m11, run: Path, release: Path, excluded: set[str]) -> dict[str, str]:
    diamond = m11.load_diamond_calls(run)
    if diamond.empty:
        return {}
    calls = m11.apply_profile_cascade(diamond, run, m11.load_cascade_config(
        release / "osmo_refdb.profile_cascade.tsv"))
    return {rid: fam for rid, fam in zip(calls["read_id"], calls["family"])
            if fam not in excluded and not fam.endswith("_decoy")}


def verify(run: Path, cat: dict[str, str | None], trans: dict[str, str]) -> tuple[int, int]:
    """(cells checked, cells that disagree) between the rebuild and gene_counts.tsv."""
    counts = next(run.glob("*.gene_counts.tsv"), None)
    if counts is None:
        return 0, 0
    rebuilt_cat = Counter(f for f in cat.values() if f)
    rebuilt_tr = Counter(trans.values())
    checked = bad = 0
    with open(counts) as fh:
        rows = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for row in rows:
            fam = row["gene"]
            pairs = [("translated_count", rebuilt_tr), ("catalog_count", rebuilt_cat)] \
                if "translated_count" in row else [("raw_count", rebuilt_tr)]
            for col, rebuilt in pairs:
                if row[col] in ("", "NA"):
                    continue
                checked += 1
                if abs(float(row[col]) - 0.5 * rebuilt.get(fam, 0)) > 1e-6:
                    bad += 1
                    print(f"  MISMATCH {run}: {fam} {col} osmotool {row[col]}, rebuilt {0.5 * rebuilt.get(fam, 0)}")
    return checked, bad


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compute-metrics", type=Path, required=True)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--annotation", type=Path, required=True)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--reads", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    m11 = load_module(args.compute_metrics)
    genes = gene_truth(args.annotation)
    excluded = load_excluded(args.release)
    args.out.mkdir(parents=True, exist_ok=True)
    summary, by_family, checks = [], [], []

    for div_dir in sorted(p for p in args.runs.iterdir() if p.is_dir()):
        div = div_dir.name
        truth = read_truth(args.reads / f"{div}.truth.tsv", genes)
        n_reads = len(truth)
        n_panel = sum(1 for f in truth.values() if f)
        for arm in ARMS:
            run = div_dir / arm
            if not run.is_dir():
                continue
            cat, decisions = catalog_calls(run) if arm != "none" else ({}, Counter())
            trans = translated_calls(m11, run, args.release, excluded)
            n_checked, n_bad = verify(run, cat, trans)
            checks.append({"divergence": div, "arm": arm, "cells_checked": n_checked, "cells_disagree": n_bad})

            c = Counter()
            fam_rows: dict[str, Counter] = defaultdict(Counter)
            for rid, fams in truth.items():
                if rid in cat:
                    call, path = cat[rid], "catalog"
                    c[f"settled_{'panel' if fams else 'nonpanel'}"] += 1
                else:
                    call, path = trans.get(rid), "translated"
                if call is None:
                    for f in fams:
                        fam_rows[f]["fn"] += 1
                    continue
                if call in fams:
                    c["correct"] += 1
                    fam_rows[call]["tp"] += 1
                else:
                    c[f"false_{path}"] += 1
                    c["false_on_nonpanel" if not fams else "false_on_other_panel"] += 1
                    fam_rows[call][f"fp_{path}"] += 1
                    for f in fams:
                        fam_rows[f]["fn"] += 1
            false = c["false_catalog"] + c["false_translated"]
            summary.append({
                "divergence": div, "arm": arm, "reads": n_reads, "panel_reads": n_panel,
                "nonpanel_reads": n_reads - n_panel,
                "settled_panel": c["settled_panel"], "settled_nonpanel": c["settled_nonpanel"],
                "ambiguous": decisions["ambiguous"],
                "correct": c["correct"], "false": false,
                "false_catalog": c["false_catalog"], "false_translated": c["false_translated"],
                "false_on_nonpanel": c["false_on_nonpanel"], "false_on_other_panel": c["false_on_other_panel"],
                "false_per_million_reads": round(false / n_reads * 1e6, 1),
                "false_per_million_nonpanel_reads": round(c["false_on_nonpanel"] / max(n_reads - n_panel, 1) * 1e6, 1),
                "recall": round(c["correct"] / n_panel, 4) if n_panel else "",
                "precision": round(c["correct"] / (c["correct"] + false), 4) if c["correct"] + false else "",
            })
            for fam in sorted(fam_rows):
                r = fam_rows[fam]
                tp, fn = r["tp"], r["fn"]
                fp = r["fp_catalog"] + r["fp_translated"]
                by_family.append({
                    "divergence": div, "arm": arm, "family": fam, "tp": tp, "fn": fn,
                    "fp_catalog": r["fp_catalog"], "fp_translated": r["fp_translated"],
                    "recall": round(tp / (tp + fn), 4) if tp + fn else "",
                    "precision": round(tp / (tp + fp), 4) if tp + fp else "",
                    "fp_per_million_reads": round(fp / n_reads * 1e6, 1),
                })

    for name, rows in (("summary.tsv", summary), ("by_family.tsv", by_family), ("verification.tsv", checks)):
        if rows:
            with open(args.out / name, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
            print(f"written: {args.out / name}")
    for r in summary:
        print(f"  d={r['divergence']:<5} {r['arm']:<5} false/M reads {r['false_per_million_reads']:>9}  "
              f"(catalog {r['false_catalog']}, translated {r['false_translated']})  recall {r['recall']}")
    if any(r["cells_disagree"] for r in checks):
        raise SystemExit("ERROR: the per-read rebuild disagrees with osmotool's gene_counts.tsv (verification.tsv)")


if __name__ == "__main__":
    main()
