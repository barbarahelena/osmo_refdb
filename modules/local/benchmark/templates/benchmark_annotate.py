#!/usr/bin/env python3
"""
benchmark_annotate.py — benchmark `osmotool annotate` on held-out
full-length proteins.

The read benchmark (BENCHMARK_SIMULATE, BENCHMARK_READS, BENCHMARK_METRICS) scores `osmotool profile`. Annotate mode
works differently: full-length ORFs, hmmsearch --cut_ga against each
family's calibrated GA cutoff (or DIAMOND blastp gated by the per-family
bitscore cutoffs from 08b), best hit per protein, decoy models competing and
then dropped. None of that is exercised by simulated reads, so it gets its
own benchmark here, at the protein level.

Evaluation set, --split test (the default; results/annotate/), all unseen
by both the models and the cutoff calibration:
  positives        every family's refs/<family>.positive.test.faa
                   (split_train_test.py: not in positive.train, which
                   builds the HMMs and DIAMOND db, and not in
                   positive.calib, which sets the cutoffs)
  hard negatives   every family's refs/<family>.negative.test.faa (same Pfam
                   domain or named paralog, different gene; negative.train
                   is what calibrated the cutoffs)
  background       refs/random_background.faa (01f: random bacterial
                   proteins without any panel Pfam domain)

--split calib (results/annotate_calib/) scores positive.calib against
negative.train instead: the calibration split, the same sequences the read
benchmark's calibration reads are simulated from. It decides what annotate
reports (--report-levels-out: osmo_refdb.annotate_report_levels.tsv), with
the report-level rule of report_rule.py at annotate's own settings
(params.annotate_*: threshold, minimum evidence, support, bootstrap), so no
report level rests on a test sequence. Each level needs bootstrap support
(protein_bootstrap: the proteins resampled, each its own source sequence;
--decision-bootstrap, --decision-support); annotate has no cascade. These are the sequences the GA and
DIAMOND cutoffs were calibrated on, so its numbers are in-sample for the
cutoffs and optimistic; the test split is the estimate to quote.

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

False positives for family F are counted as in compute_metrics.py:
  fp                  F's own hard negatives called F     -> precision, f1
  fp_cross            other families' positives called F  -> precision_panel
  fp_other_negatives  other families' hard negatives called F
  fp_all              fp + fp_cross + fp_other_negatives   -> precision_all, f1_all
  fp_background       background proteins called F (in no precision, plus
                      background_fp_per_10k)
f1_all is the report-level rule's metric, as in profile mode: it counts the
same proteins at gene and group level, where per-gene precision leaves out
the siblings' hard negatives that group precision counts.

Average precision (test split, with osmotool annotate --keep_aln's
alignments in --calls-from): per family and method, the ranking by
bit score of the proteins whose best hit is the family, over the family's
positives and its own hard negatives (report_rule.ap_row). The best hit is
taken as osmotool takes it, from every alignment it kept: hmm, hmmsearch
--cut_ga, so the ranking only reaches down to the family's GA; diamond,
every DIAMOND alignment before the per-family bit-score cutoffs, so it
reaches below them. verification.tsv checks those best hits against
osmotool's own calls.

Usage:
  python benchmark_annotate.py --refs releases/v10/refs \\
      --release releases/v10 --families families.yaml \\
      --out releases/v10/results/annotate [--methods hmm diamond] [--split calib]
Output (in --out):
  summary.tsv           per family x method
  summary_by_group.tsv  per paralog_group x method (group scored as one unit)
  overall.tsv           micro/macro totals and background rate per method
  calls.<method>.tsv    one row per evaluated protein: truth and call
  average_precision.tsv, ap_scores.<method>.tsv, verification.tsv
                        test split with kept alignments: AP per family and
                        method, the proteins it ranks, the best-hit check
"""

from __future__ import annotations

import argparse
import csv
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import yaml

import report_rule

RANDOM_BG_FASTA = "random_background.faa"
NEGATIVE_SPLIT = {"test": "test", "calib": "train"}   # positives' split -> its negatives'


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


def build_eval_set(refs: Path, families: list[str],
                   split: str = "test") -> tuple[list[dict], dict[str, int]]:
    """Assemble the evaluation proteins: positive.<split> and their hard
    negatives (negative.test, or negative.train for calib), and the
    background. Returns (rows, n_dropped_by_kind)."""
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
        pos_path = refs / f"{family}.positive.{split}.faa"
        if split == "test" and not (refs / f"{family}.positive.calib.faa").exists() and pos_path.exists():
            print(f"[{family}] WARNING: no positive.calib split -- this family's "
                  f"cutoffs were calibrated on the positives scored here",
                  file=sys.stderr)
        if not pos_path.exists():
            print(f"[{family}] WARNING: no positive.{split} split, family has no "
                  f"positives in this benchmark", file=sys.stderr)
        for header, seq in parse_fasta(pos_path):
            add("positive", family, header, seq)
        neg_path = refs / f"{family}.negative.{NEGATIVE_SPLIT[split]}.faa"
        if not neg_path.exists():
            print(f"[{family}] WARNING: no negative.{NEGATIVE_SPLIT[split]} split, family has no "
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
           "--out_prefix", str(out_prefix), "--threads", str(threads), "--keep_aln"]
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
    fp_all = fp + fp_cross + fp_other_neg
    precision_all = ratio(tp, tp + fp_all)
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
        # added after v10, as in compute_metrics.py: the report-level rule's metric
        "fp_all": fp_all, "precision_all": r3(precision_all),
        "f1_all": r3(f1_of(precision_all, recall)),
    }


# ---- the report levels (calibration split)

MANIFEST_COLUMNS = [
    "unit", "report_level", "members", "paralog_group", "method",
    "reason", "reason_detail", "rule_metric", "threshold", "n_eval_positives", "min_eval_positives",
    "tp", "fp_all", "fn", "precision_all", "recall", "f1_all", "f1_all_gene", "f1_all_group",
    # the bootstrap behind the level (protein_bootstrap); annotate has no cascade, so those two are NA
    "f1_all_p10", "f1_all_p90", "level_support", "cascade_gain", "cascade_gain_support",
]


def protein_bootstrap(rows: list[dict], calls: dict[str, str], units: dict[tuple, set[str]],
                      n_boot: int, seed: int) -> tuple[dict[tuple, list[float]], dict[tuple, tuple]]:
    """f1_all of every unit in n_boot resamples of the evaluation proteins.

    Each protein is its own source sequence, so the cluster bootstrap of the
    read benchmark (compute_metrics.cluster_bootstrap) is a protein
    bootstrap here: stratified by (source family, positive / hard
    negative), each stratum drawn with replacement to its own size, one
    draw for every unit. The background is left out, as from f1_all.
    Standard library only (this scorer's container has no numpy).

    A protein's effect on each unit (tp, fn or fp_all) is worked out once;
    a resample only counts how often each protein is drawn, and a protein
    that touches no unit (a hard negative nobody calls) is drawn but adds
    nothing. Returns ({unit key: f1_all per resample}, {unit key: (tp, fn,
    fp_all) over every protein once: the point counts, for a check against
    summary.tsv})."""
    keys = list(units)
    strata: dict[tuple, list[list[int]]] = {}
    for row in rows:
        if row["kind"] not in ("positive", "negative"):
            continue
        called = calls.get(row["protein_id"])
        effects = []
        for u, key in enumerate(keys):
            members = units[key]
            called_in = called in members
            if row["kind"] == "positive" and row["source_family"] in members:
                effects.append(3 * u + (0 if called_in else 1))
            elif called_in:
                effects.append(3 * u + 2)
        strata.setdefault((row["source_family"], row["kind"]), []).append(effects)
    plan = []
    for key in sorted(strata):
        items = strata[key]
        plan.append((len(items), [e for e in items if e]))   # proteins with an effect first: index < m

    totals = [0] * (3 * len(keys))
    for _, effectful in plan:
        for effects in effectful:
            for j in effects:
                totals[j] += 1
    point = {key: tuple(totals[3 * u:3 * u + 3]) for u, key in enumerate(keys)}

    rnd = random.Random(seed).random
    out: dict[tuple, list[float]] = {key: [] for key in keys}
    for _ in range(n_boot):
        counts = [0] * (3 * len(keys))
        for n, effectful in plan:
            m = len(effectful)
            for _ in range(n):
                i = int(rnd() * n)
                if i < m:
                    for j in effectful[i]:
                        counts[j] += 1
        for u, key in enumerate(keys):
            tp, fn, fp_all = counts[3 * u:3 * u + 3]
            den = 2 * tp + fn + fp_all
            out[key].append(2 * tp / den if den else 0.0)
    return out, point


def percentile(values: list[float], q: float) -> float:
    """numpy.percentile's default (linear interpolation), q in [0, 100]."""
    xs = sorted(values)
    pos = (len(xs) - 1) * q / 100
    lo = int(pos)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


def _value(row: dict | None, metric: str) -> float:
    """A summary row's metric as a number; no row, or no calls, counts as 0."""
    return float(row[metric]) if row and row[metric] != "" else 0.0


def decide_levels(fam_cfg: list[dict], gene_rows: dict[str, dict], group_rows: dict[str, dict],
                  method: str, min_value: float, min_n: int, metric: str = "f1_all",
                  boot: dict[tuple, list[float]] | None = None, min_support: float = 0.0) -> list[dict]:
    """One manifest row per reporting unit (a paralog group reported as a
    group, else a family), by report_rule.rule_levels on the calibration
    split's rows of one method. families.yaml has no annotate report level
    and no annotate scope, so every level here is the rule's.

    boot: {("gene", family) | ("group", group): f1_all per resample}
    (protein_bootstrap), or None: then the point estimates decide alone. A
    level needs f1_all >= min_value on the point estimate and in >=
    min_support of the resamples (report_rule.passes_level)."""
    boot = boot or {}

    def cell(row, key):
        return "NA" if not row or row[key] == "" else str(row[key])

    def counts(row):
        return {k: (row or {}).get(k, "") for k in ("tp", "fp_all", "fn", "precision_all", "recall", "f1_all")}

    def support(key):
        xs = boot.get(key)
        return None if xs is None else sum(x >= min_value for x in xs) / len(xs)

    def robust(key):
        xs = boot.get(key)
        sup = support(key)
        return {"f1_all_p10": round(percentile(xs, 10), 3) if xs else "",
                "f1_all_p90": round(percentile(xs, 90), 3) if xs else "",
                "level_support": "" if sup is None else round(sup, 3),
                "cascade_gain": "NA", "cascade_gain_support": "NA"}

    blocks, seen = [], set()
    for fam in fam_cfg:
        g = fam.get("paralog_group") or ""
        if not g:
            blocks.append(([fam["name"]], ""))
        elif g not in seen:
            seen.add(g)
            blocks.append(([f["name"] for f in fam_cfg if f.get("paralog_group") == g], g))

    rows = []
    for members, group in blocks:
        grow = group_rows.get(group) if group else None
        levels = report_rule.rule_levels(
            members, group, {f: _value(gene_rows.get(f), metric) for f in members},
            {f: (gene_rows.get(f) or {}).get("n_positive", 0) for f in members},
            _value(grow, metric) if group else None, (grow or {}).get("n_positive", 0) if group else None,
            min_value, min_n, metric,
            gene_support={f: support(("gene", f)) for f in members},
            group_support=support(("group", group)) if group else None, min_support=min_support)
        common = {"paralog_group": group, "method": method, "rule_metric": metric,
                  "threshold": min_value, "min_eval_positives": min_n,
                  "f1_all_group": cell(grow, "f1_all") if group else ""}
        if levels[members[0]][0] == "group":
            rows.append({"unit": group, "report_level": "group", "members": ",".join(members),
                         "reason": levels[members[0]][1], "reason_detail": levels[members[0]][2],
                         "n_eval_positives": grow["n_positive"], **counts(grow), **common,
                         "f1_all_gene": ",".join(cell(gene_rows.get(f), "f1_all") for f in members),
                         **robust(("group", group))})
            continue
        for f in members:
            level, reason, detail = levels[f]
            row = gene_rows.get(f)
            rows.append({"unit": f, "report_level": level, "members": f,
                         "reason": reason, "reason_detail": detail,
                         "n_eval_positives": (row or {}).get("n_positive", 0),
                         **counts(row if level == "gene" else None), **common,
                         "f1_all_gene": cell(row, "f1_all"), **robust(("gene", f))})
    return rows


def write_report_levels(path: Path, rows: list[dict], method: str, min_value: float, min_n: int,
                        n_boot: int = 0, seed: int = 0, min_support: float = 0.0) -> None:
    header = [
        f"osmo_refdb annotate report levels: one row per reporting unit, for osmotool annotate --method {method}.",
        "Decided on the calibration split (positive.calib against negative.train, full-length proteins; "
        "results/annotate_calib/), never on the test split (results/annotate/). The rule of "
        "osmo_refdb.profile_report_levels.tsv (report_rule.py), at annotate's own settings, without the cascade.",
        *report_rule.rule_header_lines(min_value, min_n, n_boot, seed, min_support, "proteins",
                                       mode="annotate"),
        "reason: meets_threshold | below_threshold | unstable | insufficient_evidence. out_of_scope and override "
        "do not occur: families.yaml sets neither a scope nor a level for annotate. reason_detail: the numbers "
        "behind it.",
        "f1_all = F1 with every non-background protein called as the unit counted as a false positive (own hard "
        "negatives, other families' positives, other families' hard negatives); see benchmark_annotate.py.",
        "tp..f1_all are calibration counts at the reported level (empty for none). f1_all_gene: per member, in "
        "members order. f1_all_group: the paralog group as one unit. cascade_gain*: NA, annotate has no cascade.",
        "The HMM GA and DIAMOND cutoffs were calibrated on these same sequences, so the numbers are in-sample "
        "for the cutoffs (optimistic); the test split is the estimate to quote.",
    ]
    with open(path, "w", newline="") as fh:
        for line in header:
            fh.write(f"# {line}\n")
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# ---- average precision (test split), from osmotool's kept alignments

ALIGNMENT_FILES = {"hmm": "annotate_hmm.hmmsearch.tblout", "diamond": "annotate_diamond.blastp.tsv"}


def read_best_hits(path: Path, method: str) -> dict[str, tuple[str, float]]:
    """{protein_id: (family, bitscore)} of each protein's best alignment in an
    osmotool annotate --keep_aln file, picked as osmotool's select_best_hits
    does (higher bitscore, then lower evalue, then smallest subject id). The
    family is the subject id up to its first '|'."""
    best: dict[str, tuple[float, float, str]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split() if method == "hmm" else line.rstrip("\n").split("\t")
            if method == "hmm":      # hmmsearch tblout: target, acc, query (the model), acc, evalue, score
                qid, sid, evalue, bits = f[0], f[2], float(f[4]), float(f[5])
            else:                    # osmotool's DIAMOND_FIELDS: qseqid sseqid pident qcovhsp scovhsp length evalue bitscore
                qid, sid, evalue, bits = f[0], f[1], float(f[6]), float(f[7])
            key = (-bits, evalue, sid)
            if qid not in best or key < best[qid]:
                best[qid] = key
    return {q: (k[2].split("|")[0], -k[0]) for q, k in best.items()}


def average_precision(rows: list[dict], hits: dict[str, tuple[str, float]], method: str,
                      families: list[str]) -> tuple[list[dict], list[dict]]:
    """Per family: AP of the best-hit bit-score ranking over its positives and
    its own hard negatives (report_rule.ap_row). A protein whose best hit is
    another family, or a decoy (osmotool drops those), is not retrieved.
    Returns (one row per family, the retrieved proteins with their scores)."""
    out, scores = [], []
    for family in families:
        pos = [r for r in rows if r["kind"] == "positive" and r["source_family"] == family]
        neg = [r for r in rows if r["kind"] == "negative" and r["source_family"] == family]
        if not pos:
            continue
        got = {"positive": [], "hard_negative": []}
        for label, items in (("positive", pos), ("hard_negative", neg)):
            for r in items:
                fam, bits = hits.get(r["protein_id"], (None, None))
                if fam == family:
                    got[label].append(bits)
                    scores.append({"method": method, "family": family, "protein_id": r["protein_id"],
                                   "accession": r["accession"], "label": label, "bitscore": bits})
        out.append({"method": method, "family": family,
                    **report_rule.ap_row(got["positive"], got["hard_negative"], len(pos), len(neg))})
    return out, scores


def verify_best_hits(rows: list[dict], hits: dict[str, tuple[str, float]],
                     calls: dict[str, str], method: str) -> dict:
    """The best hit read back from the kept alignments against osmotool's own
    call. hmm: they must agree on every protein (the alignments are the
    calls' input, decoys dropped). diamond: the per-family cutoffs are not
    applied to the best hit, so only a called protein is checked, and it
    agrees unless a cutoff dropped its top alignment."""
    n_checked = n_disagree = 0
    for r in rows:
        pid = r["protein_id"]
        fam = hits.get(pid, (None, None))[0]
        best = None if fam is None or fam.endswith("_decoy") else fam
        call = calls.get(pid)
        if method == "diamond" and call is None:
            continue
        n_checked += 1
        n_disagree += best != call
    return {"check": f"best hit from {ALIGNMENT_FILES[method]} vs gene_coordinates.tsv",
            "method": method, "proteins_checked": n_checked, "proteins_disagree": n_disagree}


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
                         "<DIR>/annotate_<method>.gene_coordinates.tsv (and, for average precision, "
                         "the alignments osmotool annotate --keep_aln kept there)")
    ap.add_argument("--split", choices=["test", "calib"], default="test",
                    help="test: positive.test against negative.test (the benchmark); calib: "
                         "positive.calib against negative.train (decides the report levels)")
    ap.add_argument("--report-levels-out", type=Path, default=None,
                    help="--split calib: write osmo_refdb.annotate_report_levels.tsv here")
    ap.add_argument("--levels-method", default="hmm", choices=["hmm", "diamond"],
                    help="the method the report levels are decided for (osmotool annotate's default: hmm)")
    ap.add_argument("--min-f1", type=float, default=0.5,
                    help="the f1_all a unit must reach on the calibration split (report_rule.py; "
                         "params.annotate_min_report_f1)")
    ap.add_argument("--decision-support", type=float, default=0.9,
                    help="the fraction of bootstrap resamples a level must hold in (params.annotate_decision_support)")
    ap.add_argument("--decision-bootstrap", type=int, default=1000,
                    help="protein-bootstrap resamples behind each level (params.annotate_decision_bootstrap; "
                         "0: point estimates alone)")
    ap.add_argument("--decision-bootstrap-seed", type=int, default=42,
                    help="the bootstrap's seed (params.decision_bootstrap_seed)")
    ap.add_argument("--min-eval-positives", type=int, default=5,
                    help="the positive calibration proteins a unit needs before the rule decides on it; "
                         "below it, none for insufficient_evidence (params.annotate_min_eval_positives). 0: off.")
    args = ap.parse_args()
    if args.report_levels_out and args.split != "calib":
        ap.error("--report-levels-out needs --split calib: report levels are never decided on the test split")
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
    rows, dropped = build_eval_set(args.refs, families, args.split)
    kinds = Counter(r["kind"] for r in rows)
    n_bg = kinds["background"]
    print(f"Evaluation set ({args.split} split): {kinds['positive']} positives "
          f"(positive.{args.split}), {kinds['negative']} hard negatives "
          f"(negative.{NEGATIVE_SPLIT[args.split]}), {n_bg} background proteins")
    for what, n in dropped.items():
        print(f"  left out: {n} {what}")
    proteins = args.out / "eval_proteins.faa"
    write_prodigal_style(rows, proteins)
    if args.prepare:
        write_dummy_assembly(args.out / "dummy_assembly.fna")
        print(f"Wrote {proteins} and {args.out / 'dummy_assembly.fna'}")
        return

    summary, group_summary, overall = [], [], []
    ap_rows, verification = [], []
    calls_by_method: dict[str, dict[str, str]] = {}
    for method in args.methods:
        print(f"\n=== osmotool annotate --method {method} ===")
        if args.calls_from:
            calls = read_calls(args.calls_from / f"annotate_{method}.gene_coordinates.tsv")
        else:
            calls = run_annotate(args.release, proteins, args.out / f"annotate_{method}",
                                 method, args.threads)
        calls_by_method[method] = calls
        unknown = sorted(set(calls.values()) - set(families))
        if unknown:
            # e.g. fused-ORF labels (mrpA_mrpB_fused): reported by osmotool,
            # not a benchmark family. Left as calls to that label, so they
            # count as misses for the true family and as nobody's FP.
            print(f"  calls to labels outside families.yaml: "
                  f"{dict(Counter(c for c in calls.values() if c in unknown))}")

        aln = (args.calls_from or args.out) / ALIGNMENT_FILES[method]
        if args.split == "test" and aln.exists():
            hits = read_best_hits(aln, method)
            verification.append(verify_best_hits(rows, hits, calls, method))
            fam_ap, scores = average_precision(rows, hits, method, families)
            ap_rows += fam_ap
            write_tsv(args.out / f"ap_scores.{method}.tsv", scores)
            v = verification[-1]
            print(f"  best hits from {aln.name}: {v['proteins_checked'] - v['proteins_disagree']}/"
                  f"{v['proteins_checked']} agree with osmotool's calls")
        elif args.split == "test":
            print(f"  no {ALIGNMENT_FILES[method]} (osmotool annotate --keep_aln): no average precision")

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
                                      "members": ",".join(sorted(members)), **s,
                                      "paralog_group": group})   # compute_metrics.py's name for it

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
    write_tsv(args.out / "average_precision.tsv", ap_rows)
    write_tsv(args.out / "verification.tsv", verification)
    print(f"\nAnnotate-mode benchmark ({args.split} split) written to {args.out}/"
          f"{{summary,summary_by_group,overall}}.tsv")

    if args.report_levels_out:
        if args.levels_method not in args.methods:
            sys.exit(f"ERROR: --levels-method {args.levels_method} was not scored (--methods)")
        gene_rows = {r["family"]: r for r in summary if r["method"] == args.levels_method}
        group_rows = {r["group"]: r for r in group_summary if r["method"] == args.levels_method}
        boot = None
        if args.decision_bootstrap > 0:
            units = {**{("gene", f): {f} for f in families},
                     **{("group", g): ms for g, ms in groups.items()}}
            t0 = time.perf_counter()
            boot, point = protein_bootstrap(rows, calls_by_method[args.levels_method], units,
                                            args.decision_bootstrap, args.decision_bootstrap_seed)
            bad = [f"{key}: summary {r['tp']}/{r['fn']}/{r['fp_all']}, bootstrap {point[key]}"
                   for key, r in ([(("gene", f), gene_rows.get(f)) for f in families]
                                  + [(("group", g), group_rows.get(g)) for g in groups])
                   if r is not None and (r["tp"], r["fn"], r["fp_all"]) != point[key]]
            if bad:
                sys.exit("ERROR: bootstrap proteins do not add up to summary.tsv:\n  " + "\n  ".join(bad))
            print(f"Decision bootstrap: {args.decision_bootstrap} resamples of the proteins, seed "
                  f"{args.decision_bootstrap_seed}, {len(units)} units, {time.perf_counter() - t0:.1f} s")
        units = decide_levels(fam_cfg, gene_rows, group_rows, args.levels_method, args.min_f1,
                              args.min_eval_positives, boot=boot, min_support=args.decision_support)
        args.report_levels_out.parent.mkdir(parents=True, exist_ok=True)
        write_report_levels(args.report_levels_out, units, args.levels_method,
                            args.min_f1, args.min_eval_positives, args.decision_bootstrap,
                            args.decision_bootstrap_seed, args.decision_support)
        print(f"Annotate report levels ({args.levels_method}) written to {args.report_levels_out}")
        for level in ("gene", "group", "none"):
            picked = [u for u in units if u["report_level"] == level]
            print(f"  {level}: {len(picked)} units: "
                  + ", ".join(f"{u['unit']}" + (f" ({u['reason']})" if level == "none" else "") for u in picked))


if __name__ == "__main__":
    main()
