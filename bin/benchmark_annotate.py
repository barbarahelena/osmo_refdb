#!/usr/bin/env python3
"""
11b_benchmark_annotate.py — benchmark `osmotool annotate` on held-out
full-length proteins.

The read benchmark (steps 9-11) scores `osmotool profile`. Annotate mode
works differently: full-length ORFs, hmmsearch --cut_ga against each
family's calibrated GA cutoff (or DIAMOND blastp gated by the per-family
bitscore cutoffs from 08b), best hit per protein, decoy models competing and
then dropped. None of that is exercised by simulated reads, so it gets its
own benchmark here, at the protein level.

Evaluation set, all unseen by both the models and the cutoff calibration:
  positives        every family's refs/<family>.positive.test.faa
                   (03_split_train_test.py: not in positive.train, which
                   builds the HMMs and DIAMOND db, and not in
                   positive.calib, which sets the cutoffs)
  hard negatives   every family's refs/<family>.negative.test.faa (same Pfam
                   domain or named paralog, different gene; negative.train
                   is what calibrated the cutoffs)
  background       refs/random_background.faa (01f: random bacterial
                   proteins without any panel Pfam domain)

A negative or background protein whose accession is a curated positive of
ANY family (train, calib or test) is left out: sibling families draw their
negatives from shared Pfam pools, so one family's negative can be another
family's true member, and calling it is not an error.

The proteins are passed to the real `osmotool annotate --proteins`, once per
method, and the per-protein calls are read back from the
<prefix>.gene_coordinates.tsv it writes. To make osmotool report
coordinates, each protein gets a Prodigal-style header and sits alone on its
own one-gene contig, so the multi-subunit system layer has nothing to
colocalise and plays no part.

False positives for family F are counted as in 11_compute_metrics.py:
  fp                  F's own hard negatives called F     -> precision, f1
  fp_cross            other families' positives called F  -> precision_panel
  fp_other_negatives  other families' hard negatives called F (reported only)
  fp_background       background proteins called F (reported only, plus
                      background_fp_per_10k)

Usage:
  python 11b_benchmark_annotate.py --refs releases/v10/refs \\
      --release releases/v10 --families families.yaml \\
      --out releases/v10/results/annotate [--methods hmm diamond]
Output (in --out):
  summary.tsv           per family x method
  summary_by_group.tsv  per paralog_group x method (group scored as one unit)
  overall.tsv           micro/macro totals and background rate per method
  calls.<method>.tsv    one row per evaluated protein: truth and call
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

RANDOM_BG_FASTA = "random_background.faa"


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records, header, seq = [], None, []
    if not path.exists():
        return records
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def accession_of(header: str) -> str:
    """'tag|accession|organism' (retag_fasta's format) -> accession; falls
    back to the first whitespace-delimited token."""
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header.split()[0]


def build_eval_set(refs: Path, families: list[str]) -> tuple[list[dict], dict[str, int]]:
    """Assemble the evaluation proteins. Returns (rows, n_dropped_by_kind)."""
    curated_positive: set[str] = set()
    for family in families:
        for header, _ in parse_fasta(refs / f"{family}.positive.faa"):
            curated_positive.add(accession_of(header))

    rows: list[dict] = []
    dropped: Counter = Counter()

    def add(kind: str, source: str, header: str, seq: str) -> None:
        rows.append({"kind": kind, "source_family": source,
                     "accession": accession_of(header), "seq": seq.rstrip("*")})

    for family in families:
        pos_path = refs / f"{family}.positive.test.faa"
        if not (refs / f"{family}.positive.calib.faa").exists() and pos_path.exists():
            print(f"[{family}] WARNING: no positive.calib split -- this family's "
                  f"cutoffs were calibrated on the positives scored here",
                  file=sys.stderr)
        for header, seq in parse_fasta(pos_path):
            add("positive", family, header, seq)
        neg_path = refs / f"{family}.negative.test.faa"
        if not neg_path.exists():
            print(f"[{family}] WARNING: no negative.test split, family has no "
                  f"hard negatives in this benchmark", file=sys.stderr)
        for header, seq in parse_fasta(neg_path):
            if accession_of(header) in curated_positive:
                dropped["hard negatives that are a curated positive of some family"] += 1
                continue
            add("negative", family, header, seq)
    for header, seq in parse_fasta(refs / RANDOM_BG_FASTA):
        if accession_of(header) in curated_positive:
            dropped["background proteins that are a curated positive"] += 1
            continue
        add("background", "random_bg", header, seq)

    for i, row in enumerate(rows):
        row["contig"] = f"s{i:07d}"
        row["protein_id"] = f"s{i:07d}_1"
    return rows, dict(dropped)


def write_prodigal_style(rows: list[dict], path: Path) -> None:
    with open(path, "w") as fh:
        for row in rows:
            end = 3 * len(row["seq"]) + 3
            fh.write(f">{row['protein_id']} # 1 # {end} # 1 # "
                     f"ID={row['contig']}_1;partial=00;start_type=ATG\n")
            seq = row["seq"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def run_annotate(release: Path, proteins: Path, out_prefix: Path, method: str,
                 threads: int) -> dict[str, str]:
    """Run osmotool annotate; return {protein_id: called family}."""
    dummy = out_prefix.parent / "dummy_assembly.fna"
    write_dummy_assembly(dummy)
    cmd = ["osmotool", "annotate", str(release), str(dummy),
           "--proteins", str(proteins), "--method", method,
           "--out_prefix", str(out_prefix), "--threads", str(threads)]
    log_path = Path(str(out_prefix) + ".log")
    with open(log_path, "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    coords = Path(str(out_prefix) + ".gene_coordinates.tsv")
    if not coords.exists():
        sys.exit(f"ERROR: {coords} not written -- see {log_path}")
    return read_calls(coords)


def write_dummy_assembly(path: Path) -> None:
    """annotate needs an assembly argument even with --proteins."""
    path.write_text(">dummy\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")


def read_calls(coords: Path) -> dict[str, str]:
    """{protein_id: called family} from osmotool annotate's gene_coordinates.tsv."""
    if not coords.exists():
        sys.exit(f"ERROR: {coords} not found")
    calls: dict[str, str] = {}
    with open(coords) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["protein_id"] in calls:
                sys.exit(f"ERROR: {row['protein_id']} called twice in {coords}")
            calls[row["protein_id"]] = row["gene"]
    return calls


def ratio(num: float, den: float):
    return num / den if den else float("nan")


def f1_of(p: float, r: float):
    if p != p or r != r or not (p + r):
        return float("nan")
    return 2 * p * r / (p + r)


def r3(x):
    return round(x, 3) if x == x else ""


def score(rows: list[dict], calls: dict[str, str], members: set[str],
          n_background: int) -> dict:
    """Score one family, or one paralog group as a unit (members = its
    families): a call to any member counts for the group."""
    tp = fn = fp = fp_cross = fp_other_neg = fp_bg = n_neg = swaps = 0
    for row in rows:
        called = calls.get(row["protein_id"])
        in_group = row["source_family"] in members
        called_in = called in members
        if row["kind"] == "positive":
            if in_group:
                if called_in:
                    tp += 1
                    swaps += called != row["source_family"]
                else:
                    fn += 1
            elif called_in:
                fp_cross += 1
        elif row["kind"] == "negative":
            if in_group:
                n_neg += 1
                fp += called_in
            elif called_in:
                fp_other_neg += 1
        elif called_in:
            fp_bg += 1
    precision = ratio(tp, tp + fp)
    precision_panel = ratio(tp, tp + fp + fp_cross)
    recall = ratio(tp, tp + fn)
    return {
        "n_positive": tp + fn, "n_hard_negative": n_neg,
        "tp": tp, "fn": fn, "fp": fp, "fp_cross": fp_cross,
        "fp_other_negatives": fp_other_neg, "fp_background": fp_bg,
        "precision": r3(precision), "recall": r3(recall),
        "f1": r3(f1_of(precision, recall)),
        "precision_panel": r3(precision_panel),
        "f1_panel": r3(f1_of(precision_panel, recall)),
        "hard_negative_fp_rate": r3(ratio(fp, n_neg)),
        "background_fp_per_10k": (round(fp_bg / n_background * 1e4, 2)
                                  if n_background else ""),
        "within_group_swaps": swaps,
    }


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", required=True, type=Path)
    ap.add_argument("--release", type=Path, default=None)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--methods", nargs="+", default=["hmm", "diamond"],
                    choices=["hmm", "diamond"])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--prepare", action="store_true",
                    help="only write <out>/eval_proteins.faa and <out>/dummy_assembly.fna, for "
                         "running osmotool annotate elsewhere (the Nextflow pipeline does)")
    ap.add_argument("--calls-from", type=Path, default=None,
                    help="do not run osmotool: read the calls from "
                         "<DIR>/annotate_<method>.gene_coordinates.tsv")
    args = ap.parse_args()
    if args.release is None and not (args.prepare or args.calls_from):
        ap.error("--release is required unless --prepare or --calls-from is given")

    with open(args.families) as fh:
        fam_cfg = yaml.safe_load(fh)["families"]
    families = [f["name"] for f in fam_cfg]
    groups: dict[str, set[str]] = defaultdict(set)
    for f in fam_cfg:
        if f.get("paralog_group"):
            groups[f["paralog_group"]].add(f["name"])

    args.out.mkdir(parents=True, exist_ok=True)
    rows, dropped = build_eval_set(args.refs, families)
    kinds = Counter(r["kind"] for r in rows)
    n_bg = kinds["background"]
    print(f"Evaluation set: {kinds['positive']} held-out positives, "
          f"{kinds['negative']} held-out hard negatives, {n_bg} background proteins")
    for what, n in dropped.items():
        print(f"  left out: {n} {what}")
    proteins = args.out / "eval_proteins.faa"
    write_prodigal_style(rows, proteins)
    if args.prepare:
        write_dummy_assembly(args.out / "dummy_assembly.fna")
        print(f"Wrote {proteins} and {args.out / 'dummy_assembly.fna'}")
        return

    summary, group_summary, overall = [], [], []
    for method in args.methods:
        print(f"\n=== osmotool annotate --method {method} ===")
        if args.calls_from:
            calls = read_calls(args.calls_from / f"annotate_{method}.gene_coordinates.tsv")
        else:
            calls = run_annotate(args.release, proteins, args.out / f"annotate_{method}",
                                 method, args.threads)
        unknown = sorted(set(calls.values()) - set(families))
        if unknown:
            # e.g. fused-ORF labels (mrpA_mrpB_fused): reported by osmotool,
            # not a benchmark family. Left as calls to that label, so they
            # count as misses for the true family and as nobody's FP.
            print(f"  calls to labels outside families.yaml: "
                  f"{dict(Counter(c for c in calls.values() if c in unknown))}")

        write_tsv(args.out / f"calls.{method}.tsv", [
            {"protein_id": r["protein_id"], "kind": r["kind"],
             "source_family": r["source_family"], "accession": r["accession"],
             "called": calls.get(r["protein_id"], "NO_CALL")} for r in rows])

        fam_rows = []
        for family in families:
            s = score(rows, calls, {family}, n_bg)
            if s["n_positive"] == 0:
                continue
            s.pop("within_group_swaps")
            fam_rows.append({"method": method, "family": family, **s})
        summary += fam_rows
        for group, members in sorted(groups.items()):
            s = score(rows, calls, members, n_bg)
            if s["n_positive"]:
                group_summary.append({"method": method, "group": group,
                                      "members": ",".join(sorted(members)), **s})

        tp = sum(r["tp"] for r in fam_rows)
        fn = sum(r["fn"] for r in fam_rows)
        fp = sum(r["fp"] for r in fam_rows)
        fp_cross = sum(r["fp_cross"] for r in fam_rows)
        n_bg_called = sum(1 for r in rows if r["kind"] == "background"
                          and calls.get(r["protein_id"]) in families)
        n_neg_called = sum(1 for r in rows if r["kind"] == "negative"
                           and calls.get(r["protein_id"]) in families)

        def mean(col):
            vals = [r[col] for r in fam_rows if r[col] != ""]
            return r3(sum(vals) / len(vals)) if vals else ""

        def median(col):
            vals = sorted(r[col] for r in fam_rows if r[col] != "")
            if not vals:
                return ""
            mid = len(vals) // 2
            return r3(vals[mid] if len(vals) % 2 else (vals[mid - 1] + vals[mid]) / 2)

        micro_p, micro_r = ratio(tp, tp + fp), ratio(tp, tp + fn)
        overall.append({
            "method": method, "n_families": len(fam_rows),
            "n_positive": tp + fn, "n_hard_negative": kinds["negative"],
            "n_background": n_bg,
            "micro_precision": r3(micro_p), "micro_recall": r3(micro_r),
            "micro_f1": r3(f1_of(micro_p, micro_r)),
            "micro_precision_panel": r3(ratio(tp, tp + fp + fp_cross)),
            "macro_precision": mean("precision"), "macro_recall": mean("recall"),
            "macro_f1": mean("f1"), "median_f1": median("f1"),
            "median_recall": median("recall"),
            "hard_negatives_called_any_family": n_neg_called,
            "hard_negative_any_call_rate": r3(ratio(n_neg_called, kinds["negative"])),
            "background_called_any_family": n_bg_called,
            "background_fp_per_10k": round(n_bg_called / n_bg * 1e4, 2) if n_bg else "",
        })
        print(f"  micro precision {r3(micro_p)}  recall {r3(micro_r)}  "
              f"median F1 {median('f1')}  background calls {n_bg_called}/{n_bg}")

    write_tsv(args.out / "summary.tsv", summary)
    write_tsv(args.out / "summary_by_group.tsv", group_summary)
    write_tsv(args.out / "overall.tsv", overall)
    print(f"\nAnnotate-mode benchmark written to {args.out}/"
          f"{{summary,summary_by_group,overall}}.tsv")


if __name__ == "__main__":
    main()
