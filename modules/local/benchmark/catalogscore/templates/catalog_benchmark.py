#!/usr/bin/env python3
"""
catalog_benchmark.py -- score `osmotool profile --catalog` on the release's
simulated test reads, with two catalogs.

A catalog changes which reads reach translated search, so its value depends
on whether a sample's strains are in the catalog. Two catalogs are scored on
the same test reads as the read benchmark:

  A  catalogued      the catalog holds every sequence of the release: train,
                     calib and test positives, both negative splits and the
                     random background. Every test read's own gene is in it:
                     the upper bound, where the catalog settles what it can.
  B  not catalogued  the catalog holds only what the benchmark does not score
                     (positive.train, positive.calib, negative.train). No test
                     read's own gene is in it, so no read has an exact gene
                     match. B is not a lower bound: it still settles reads
                     through conserved regions of homologous genes.

Why B settles reads. The catalog's identity is 1 - NM / aligned columns: the
identity of the stretch the read covers, not of the whole gene. CD-HIT at 90%
keeps every catalog gene at most about 90% identical, end to end, to a test
gene, but two such genes still contain conserved windows that are >= 97%
identical, and a read that lands in one passes the 97% bar. The reads B
settles are mostly called as the right family (settled_reads.tsv): B assigns
reads from strains whose own gene is absent from the catalog. A catalog with
no homologous window at all would be a true lower bound; for this panel it probably
does not exist. The reference for "no catalog" is the diamond_cascade line.

Both catalogs are indexed as the release's catalog is (build-catalog --index
neighbours: the genes with an HMM hit and the genes close to them).

The catalog genes are the release's proteins back-translated exactly as the
read constructs are (BENCHMARK_CATALOGGENES), so a test read matches its own
gene base for base. That same choice makes simulated nucleotide identity
equal protein identity, which real strains do not match: they also differ at
synonymous sites, so their windows reach 97% less often. This applies to B
too, whose settled reads are therefore more than real uncatalogued strains
would give. It is also why only these two catalogs are scored, not a "related
strain" scenario.

Per read, the hybrid call is the catalog's decision when the catalog settles
the read (a unique call counts, a drop or an ambiguous read is no call) and
otherwise the read's translated call (DIAMOND + the release's cascade, from
the test benchmark). That equals what profile --catalog reports, because
every step decides each read on its own; BENCHMARK_CATALOGMAP's full
profile --catalog runs on a few samples check it (verification.tsv).
Methods scored (summary_by_report_unit.tsv, same metrics and units as the
read benchmark):

  diamond_cascade          no catalog (the read benchmark's numbers, recomputed)
  catalog_<S>              hybrid with scenario S
  catalog_only_<S>         catalog calls only, no translated fallback
                           (--catalog_only): what a catalog-only profiler sees

mixture.tsv: f1_all per reporting unit for samples where a fraction f of the
reads comes from catalogued strains (counts mixed linearly between A and B).
It interpolates between the two scored catalogs; neither end is a bound
(A holds each read's own gene, where real catalogs hold relatives; B settles
reads through conserved windows, see above).

Usage:
  python catalog_benchmark.py --metrics-module compute_metrics.py --reads reads \\
      --diamond run/diamond --cascade cascade --families families.yaml \\
      --runs catalog_runs --out catalog_benchmark
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import statistics as st
from pathlib import Path

import pandas as pd

SCENARIOS = ["A", "B"]
MIXTURE_F = [0.0, 0.25, 0.5, 0.75, 1.0]


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("compute_metrics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def samples(reads: Path) -> list[str]:
    return sorted(p.name[:-len("_R1.fastq")] for p in reads.glob("*_R1.fastq"))


def catalog_decisions(run_dir: Path, sample_names: list[str]) -> pd.DataFrame:
    """Every settled read of every sample: read_id (with /1, /2 as in
    read_truth.tsv), sample, decision, label."""
    frames = []
    for s in sample_names:
        f = run_dir / f"{s}.catalog_hits.tsv"
        if not f.exists() or f.stat().st_size == 0:
            continue
        df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        if df.empty:
            continue
        df["read_id"] = df["read"] + "/" + df["mate"]
        df["sample"] = s
        frames.append(df[["read_id", "sample", "decision", "label"]])
    if not frames:
        return pd.DataFrame(columns=["read_id", "sample", "decision", "label"])
    return pd.concat(frames, ignore_index=True)


def hybrid_calls(translated: pd.DataFrame, settled: pd.DataFrame, *, fallback: bool) -> pd.DataFrame:
    """Per-read calls of profile --catalog: catalog calls for settled reads,
    translated calls for the rest (none if not fallback)."""
    calls = settled[settled["decision"] == "call"].rename(columns={"label": "family"})
    calls = calls[["read_id", "sample", "family"]].assign(bitscore=1e9, evalue=0.0)
    if not fallback:
        return calls.reset_index(drop=True)
    key = set(zip(settled["sample"], settled["read_id"]))
    keep = [k not in key for k in zip(translated["sample"], translated["read_id"])]
    rest = translated[keep][["read_id", "sample", "family", "bitscore", "evalue"]]
    return pd.concat([rest, calls], ignore_index=True)


def verify(verify_dir: Path, translated: pd.DataFrame, settled: pd.DataFrame) -> list[dict]:
    """The full profile --catalog runs (with fallback, cascade on) versus
    the per-read rebuild."""
    rows = []
    for counts in sorted(verify_dir.glob("*.gene_counts.tsv")):
        s = counts.name[:-len(".gene_counts.tsv")]
        reported = pd.read_csv(counts, sep="\t", comment="#")
        mine = hybrid_calls(translated[translated["sample"] == s], settled[settled["sample"] == s],
                            fallback=True)
        rebuilt = (mine.groupby("family").size() * 0.5).to_dict()
        n_bad = 0
        for fam, raw in zip(reported["gene"], reported["raw_count"]):
            if abs(rebuilt.get(fam, 0.0) - raw) > 1e-6:
                n_bad += 1
                print(f"  MISMATCH {s} {fam}: osmotool {raw}, rebuilt {rebuilt.get(fam, 0.0)}")
        rows.append({"sample": s, "families": len(reported), "disagree": n_bad})
    return rows


def mixture(m11, rows_a: dict, rows_b: dict, units: list[dict]) -> list[dict]:
    out = []
    for u in units:
        a, b = rows_a.get(u["unit"]), rows_b.get(u["unit"])
        if a is None or b is None:
            continue
        row = {"unit": u["unit"], "report_level": u["report_level"]}
        for f in MIXTURE_F:
            tp = f * a["tp"] + (1 - f) * b["tp"]
            fp = f * a["fp_all"] + (1 - f) * b["fp_all"]
            fn = f * a["fn"] + (1 - f) * b["fn"]
            p = tp / (tp + fp) if tp + fp else float("nan")
            r = tp / (tp + fn) if tp + fn else float("nan")
            row[f"f1_all_f{f:g}"] = round(m11._f1(p, r), 3) if p == p and r == r else ""
        out.append(row)
    return out


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Written: {path}")


def counts_of(rows: dict) -> dict:
    return {u: {k: float(v) if v != "" else 0.0 for k, v in r.items() if k in ("tp", "fp_all", "fn")}
            for u, r in rows.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics-module", type=Path, required=True, help="compute_metrics.py")
    ap.add_argument("--reads", type=Path, required=True, help="the test reads, with read_truth.tsv")
    ap.add_argument("--diamond", type=Path, required=True, help="the test benchmark's DIAMOND results")
    ap.add_argument("--cascade", type=Path, required=True,
                    help="directory with osmo_refdb.profile_cascade.tsv and ..._report_levels.tsv")
    ap.add_argument("--families", type=Path, required=True)
    ap.add_argument("--runs", type=Path, required=True,
                    help="one <scenario>/ per catalog: runs/<sample>.catalog_hits.tsv, verify/")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    m11 = load_module(args.metrics_module)
    args.out.mkdir(parents=True, exist_ok=True)
    meta = m11.load_family_meta(args.families)
    families = list(meta)
    groups = {f: m["group"] for f, m in meta.items() if m["group"]}
    units = m11.load_report_levels(args.cascade / "osmo_refdb.profile_report_levels.tsv")
    names = samples(args.reads)

    truth = m11.load_read_truth(args.reads)
    diamond = m11.load_diamond_calls(args.diamond)
    translated = m11.apply_profile_cascade(
        diamond, args.diamond, m11.load_cascade_config(args.cascade / "osmo_refdb.profile_cascade.tsv"))
    unit_rows, family_rows, bg_rows, settled_rows = [], [], [], []
    by_method: dict[str, dict] = {}

    def score(calls: pd.DataFrame, method: str) -> None:
        merged = m11.label_reads(calls, truth)
        rows = m11.summarize_report_units(merged, method, units)
        unit_rows.extend(rows)
        by_method[method] = {r["unit"]: r for r in rows}
        family_rows.extend(m11.summarize(merged, method, families, groups))
        bg_rows.extend(m11.summarize_background_fpr(merged, method)[0])

    score(translated, "diamond_cascade")
    settled_by = {}
    t = truth.set_index("read_id")
    for name in SCENARIOS:
        settled = catalog_decisions(args.runs / name / "runs", names)
        settled_by[name] = settled
        score(hybrid_calls(translated, settled, fallback=True), f"catalog_{name}")
        if name != "B":
            score(hybrid_calls(translated, settled, fallback=False), f"catalog_only_{name}")
        lab = settled.join(t[["family", "label"]], on="read_id", rsuffix="_true")
        for (fam_true, label), grp in lab.groupby(["family", "label_true"]):
            settled_rows.append({
                "scenario": name, "source": fam_true, "kind": label, "n_settled": len(grp),
                "call": int((grp["decision"] == "call").sum()),
                "call_correct_family": int(((grp["decision"] == "call") &
                                            (grp["label"] == fam_true)).sum()),
                "drop": int((grp["decision"] == "drop").sum()),
                "ambiguous": int((grp["decision"] == "ambiguous").sum()),
            })

    mix = mixture(m11, counts_of(by_method["catalog_A"]), counts_of(by_method["catalog_B"]), units)

    print("Verifying the per-read rebuild against the full profile --catalog runs", flush=True)
    ver = verify(args.runs / "A" / "verify", translated, settled_by["A"])

    write_tsv(args.out / "summary_by_report_unit.tsv", unit_rows)
    write_tsv(args.out / "summary.tsv", family_rows)
    write_tsv(args.out / "background_fpr.tsv", bg_rows)
    write_tsv(args.out / "settled_reads.tsv", settled_rows)
    write_tsv(args.out / "mixture.tsv", mix)
    write_tsv(args.out / "verification.tsv", ver)

    reported = [u["unit"] for u in units if u["report_level"] != "none"]

    def med(method, key="f1_all"):
        vals = [float(by_method[method][u][key]) if by_method[method][u][key] != "" else 0.0
                for u in reported if u in by_method[method]]
        return st.median(vals) if vals else float("nan")
    print("\nMedian over the reported units (f1_all / precision_all / recall):")
    for method in by_method:
        print(f"  {method:22} {med(method):.3f} / {med(method, 'precision_all'):.3f} / "
              f"{med(method, 'recall'):.3f}")
    for name, settled in settled_by.items():
        print(f"  {name}: {len(settled)} test reads settled by the catalog")
    for r in ver:
        print(f"  verify {r['sample']}: {r['families'] - r['disagree']}/{r['families']} families agree")
    if any(r["disagree"] for r in ver):
        raise SystemExit("ERROR: the per-read rebuild disagrees with profile --catalog (verification.tsv)")


if __name__ == "__main__":
    main()
