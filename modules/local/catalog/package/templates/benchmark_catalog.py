#!/usr/bin/env python3
"""
benchmark_catalog.py -- benchmark `osmotool profile --catalog` on the
release's simulated test reads, in two bounding scenarios.

A catalog changes which reads reach translated search, so its value depends
on whether a sample's strains are in the catalog. Two bounds are scored on
the same test reads as steps 9-11:

  A  catalogued      the catalog holds every sequence of the release: train,
                     calib and test positives, both negative splits and the
                     random background. Every test read's own gene is in it:
                     the upper bound, where the catalog settles what it can.
  B  not catalogued  the catalog holds only what the benchmark does not score
                     (positive.train, positive.calib, negative.train). No test
                     read's own gene is in it; CD-HIT at 90% keeps every
                     catalog gene <= ~90% identical to a test gene, below the
                     97% bar, so B should settle (almost) nothing and
                     reproduce the translated-search numbers. Any test read B
                     settles anyway is reported.

Both are indexed in full (index all, competitive). A is also scored with the
small index a release ships (A_hits: genes with an HMM hit only).

Catalog genes are the release's proteins back-translated exactly as
simulate_reads.py builds its read constructs (protein_to_dna: one fixed
codon per amino acid), so a test read matches its own gene base for base.
That same choice makes simulated nucleotide identity track protein identity,
which real strains do not: in B it is harmless (no catalog gene reaches 97%),
but it is why only these two bounds are scored, not a "related strain"
scenario.

Per read, the hybrid call is the catalog's decision when the catalog settles
the read (a unique call counts, a drop or an ambiguous read is no call) and
otherwise the read's call from steps 10-11 (DIAMOND + the release's cascade).
That equals what profile --catalog reports, because every step decides each
read on its own; this is checked by running the full profile --catalog on a
few samples (verification.tsv). Methods scored (summary_by_report_unit.tsv,
same metrics and units as step 11):

  diamond_cascade          no catalog (step 11's numbers, recomputed)
  catalog_<S>              hybrid with scenario S
  catalog_only_<S>         catalog calls only, no translated fallback
                           (--catalog_only): what a catalog-only profiler sees

mixture.tsv: f1_all per reporting unit for samples where a fraction f of the
reads comes from catalogued strains (counts mixed linearly between A_full and
B_full).

Usage:
  python bin/benchmark_catalog.py --release releases/v10 \\
      --families families.yaml --threads 8
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import statistics as st
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)   # an absolute path stays itself
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# the benchmark's own scripts, in their modules (modules/local/benchmark/<step>/templates/)
BENCHMARK = HERE.parents[2] / "benchmark"
m11 = _load("compute_metrics", BENCHMARK / "metrics" / "templates" / "compute_metrics.py")
m09 = _load("simulate_reads", BENCHMARK / "simulate" / "templates" / "simulate_reads.py")

SCENARIO_SPLITS = {
    "A": ["positive.train", "positive.calib", "positive.test", "negative.train", "negative.test"],
    "B": ["positive.train", "positive.calib", "negative.train"],
}
SCENARIOS = {"A_full": ("A", "all"), "A_hits": ("A", "hits"), "B_full": ("B", "all")}
VERIFY_SAMPLES = ["opuAC.negative", "proX.positive", "trkH.negative"]
MIXTURE_F = [0.0, 0.25, 0.5, 0.75, 1.0]


def parse_fasta(path: Path):
    return m09.parse_fasta(path)


def write_genes(refs: Path, families: list[str], scenario: str, out: Path) -> int:
    """Back-translated catalog genes for a scenario; IDs carry the source
    file so a protein present in two files stays two genes."""
    n = 0
    with open(out, "w") as fh:
        files = [refs / f"{fam}.{split}.faa" for fam in families
                 for split in SCENARIO_SPLITS[scenario]]
        if scenario == "A":
            files.append(refs / m09.RANDOM_BG_FASTA)
        for f in files:
            if not f.exists():
                continue
            tag = f.name[:-4]
            for header, protein in parse_fasta(f):
                fh.write(f">{tag}|{header.split()[0]}\n{m09.protein_to_dna(protein)}\n")
                n += 1
    return n


def check_constructs(reads: Path, refs: Path, families: list[str]) -> int:
    """The test constructs' CDS (between the 300 bp flanks) must be
    protein_to_dna of the test proteins, or the catalog genes would not be
    the read sources."""
    checked = 0
    for fam in families[:5]:
        ctx = dict(parse_fasta(reads / "genomic_context" / f"{fam}.positive.fasta"))
        for header, protein in parse_fasta(refs / f"{fam}.positive.test.faa")[:20]:
            cds = ctx[header][m09.FLANK_LEN:len(ctx[header]) - m09.FLANK_LEN]
            if cds != m09.protein_to_dna(protein):
                sys.exit(f"ERROR: construct of {header} is not protein_to_dna of its protein")
            checked += 1
    return checked


def run(cmd: list[str], log: Path) -> None:
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        if subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode != 0:
            sys.exit(f"ERROR: command failed, see {log}: {' '.join(cmd)}")


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


def verify(release: Path, reads: Path, catalog: Path, out: Path, translated: pd.DataFrame,
           settled: pd.DataFrame, threads: int, log: Path) -> list[dict]:
    """Full profile --catalog (with fallback, cascade on) on a few samples
    versus the per-read rebuild."""
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in VERIFY_SAMPLES:
        r1, r2 = reads / f"{s}_R1.fastq", reads / f"{s}_R2.fastq"
        if not r1.exists():
            continue
        run(["osmotool", "profile", str(release), "-1", str(r1), "-2", str(r2),
             "--out_prefix", str(out / s), "--catalog", str(catalog), "--threads", str(threads),
             "--min_identity", "0.40", "--min_query_cover", "0.50", "--no_markers"], log)
        reported = pd.read_csv(out / f"{s}.gene_counts.tsv", sep="\t", comment="#")
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


def mixture(rows_a: dict, rows_b: dict, units: list[dict]) -> list[dict]:
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None,
                    help="default: <release>/results/catalog_benchmark")
    args = ap.parse_args()

    rel, refs = args.release, args.release / "refs"
    results = rel / "results"
    reads = results / "reads"
    out = args.out or results / "catalog_benchmark"
    out.mkdir(parents=True, exist_ok=True)
    log = out / "commands.log"
    meta = m11.load_family_meta(args.families)
    families = list(meta)
    groups = {f: m["group"] for f, m in meta.items() if m["group"]}
    units = m11.load_report_levels(rel / "osmo_refdb.profile_report_levels.tsv")
    names = samples(reads)

    print(f"Constructs checked against protein_to_dna: {check_constructs(reads, refs, families)}")

    # --- catalogs ---
    for scen in ("A", "B"):
        genes = out / f"genes_{scen}.fna"
        if not genes.exists():
            print(f"Scenario {scen}: {write_genes(refs, families, scen, genes)} catalog genes")
    for name, (scen, index) in SCENARIOS.items():
        cdir = out / "catalogs" / name
        if not (cdir / "catalog.json").exists():
            print(f"Building catalog {name} (index {index})", flush=True)
            run(["osmotool", "build-catalog", str(rel), "--genes", str(out / f"genes_{scen}.fna"),
                 "-o", str(cdir), "--index", index, "-t", str(args.threads)], log)

    # --- catalog decisions on every test sample (catalog only: no DIAMOND) ---
    for name in SCENARIOS:
        rdir = out / "runs" / name
        rdir.mkdir(parents=True, exist_ok=True)
        todo = [s for s in names if not (rdir / f"{s}.catalog_hits.tsv").exists()]
        print(f"Mapping {len(todo)} samples to {name}", flush=True)
        for s in todo:
            run(["osmotool", "profile", str(rel), "-1", str(reads / f"{s}_R1.fastq"),
                 "-2", str(reads / f"{s}_R2.fastq"), "--out_prefix", str(rdir / s),
                 "--catalog", str(out / "catalogs" / name), "--catalog_only", "--no_cascade",
                 "--no_markers", "--keep_aln", "--threads", str(args.threads)], log)

    # --- scoring ---
    truth = m11.load_read_truth(reads)
    diamond = m11.load_diamond_calls(results / "diamond")
    translated = m11.apply_profile_cascade(
        diamond, results / "diamond", m11.load_cascade_config(rel / "osmo_refdb.profile_cascade.tsv"))
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
    for name in SCENARIOS:
        settled = catalog_decisions(out / "runs" / name, names)
        settled_by[name] = settled
        score(hybrid_calls(translated, settled, fallback=True), f"catalog_{name}")
        if name != "B_full":
            score(hybrid_calls(translated, settled, fallback=False), f"catalog_only_{name}")
        t = truth.set_index("read_id")
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

    mix = mixture({u: {k: float(v) if v != "" else 0.0 for k, v in r.items() if k in ("tp", "fp_all", "fn")}
                   for u, r in by_method["catalog_A_full"].items()},
                  {u: {k: float(v) if v != "" else 0.0 for k, v in r.items() if k in ("tp", "fp_all", "fn")}
                   for u, r in by_method["catalog_B_full"].items()}, units)

    print("Verifying the per-read rebuild against full profile --catalog runs", flush=True)
    ver = verify(rel, reads, out / "catalogs" / "A_full", out / "verify", translated,
                 settled_by["A_full"], args.threads, log)

    mdir = out / "metrics"
    mdir.mkdir(exist_ok=True)
    write_tsv(mdir / "summary_by_report_unit.tsv", unit_rows)
    write_tsv(mdir / "summary.tsv", family_rows)
    write_tsv(mdir / "background_fpr.tsv", bg_rows)
    write_tsv(mdir / "settled_reads.tsv", settled_rows)
    write_tsv(mdir / "mixture.tsv", mix)
    write_tsv(mdir / "verification.tsv", ver)

    reported = [u["unit"] for u in units if u["report_level"] != "none"]

    def med(method, key="f1_all", which=reported):
        vals = [float(by_method[method][u][key]) if by_method[method][u][key] != "" else 0.0
                for u in which if u in by_method[method]]
        return st.median(vals) if vals else float("nan")
    print("\nMedian over the reported units (f1_all / precision_all / recall):")
    for method in by_method:
        print(f"  {method:22} {med(method):.3f} / {med(method, 'precision_all'):.3f} / "
              f"{med(method, 'recall'):.3f}")
    for name, settled in settled_by.items():
        print(f"  {name}: {len(settled)} test reads settled by the catalog")
    for r in ver:
        print(f"  verify {r['sample']}: {r['families'] - r['disagree']}/{r['families']} families agree")


if __name__ == "__main__":
    main()
