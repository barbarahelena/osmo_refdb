#!/usr/bin/env python3
"""
11_compute_metrics.py — score DIAMOND vs. HMM against the simulated-read truth table

Parses:
  - results/diamond/<sample>.blastx.tsv (osmotool --keep_aln output)
  - results/hmm/<sample>.hmmscan.tblout (hmmscan output, raw bit scores)
  - results/reads/read_truth.tsv (from 09b_compute_read_truth.py; per-read
    truth based on actual CDS overlap, not just source-construct label)

For each family and method, computes read-level precision/recall/F1. A read
is reduced to its single best-scoring family call (see best_hit_per_read) so
that one read can only ever count as one false positive, not once per
spurious model hit -- important for promiscuous domains (e.g. ectB, betL)
that legitimately score weakly against many unrelated HMMs/DIAMOND targets.

Also sweeps DIAMOND's/HMM's bit-score thresholds to produce ROC/PR-style
curves (saved as PNGs) plus two summary tables:
  - summary.tsv: precision/recall/F1 accepting every best-hit call as-is
  - best_threshold_summary.tsv: precision/recall/F1 at the per-family score
    threshold that maximizes F1 (i.e. the best achievable operating point if
    a stricter score cutoff were applied on top of the best-hit call)

Usage:
  python 11_compute_metrics.py --results results/ --out results/metrics --families families.yaml
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from sklearn.metrics import precision_recall_curve, roc_curve, auc


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def load_read_truth(reads_dir: Path) -> pd.DataFrame:
    """
    Load per-read ground truth from 09b_compute_read_truth.py's read_truth.tsv
    (read_id, family, label, uniprot_id, overlaps_cds, overlap_bp). Unlike the
    contig-level truth.tsv (one row per source protein), this reflects whether
    each individual simulated read actually overlaps the true CDS -- a read
    from a "positive" construct that lands entirely in the random flanking
    DNA carries no gene signal and must not be scored as a true positive.
    """
    return pd.read_csv(reads_dir / "read_truth.tsv", sep="\t")


def load_diamond_calls(diamond_dir: Path) -> pd.DataFrame:
    """
    Parse osmotool's kept alignment TSVs (--keep_aln) to get per-read best-hit
    family assignments + score (bitscore) for threshold sweeping.
    """
    rows = []
    for aln_file in diamond_dir.glob("*.blastx.tsv"):
        if aln_file.stat().st_size == 0:
            continue
        sample = aln_file.stem.replace(".blastx", "")
        df = pd.read_csv(aln_file, sep="\t", header=None)
        # osmotool's custom diamond --outfmt: see DIAMOND_FIELDS in
        # src/osmotool/functions.py (NOT the 12-column BLAST default).
        # Includes scovhsp (subject coverage) since the subject-coverage
        # filter was added -- keep this column list in sync with that file.
        df.columns = [
            "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "length", "evalue", "bitscore",
        ][: df.shape[1]]
        # keep only the best (highest-bitscore) hit per read for scoring
        df = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
        df["sample"] = sample
        df["family"] = df["sseqid"].str.split("|").str[0]
        # osmotool reports qseqid as the raw read ID, unchanged
        df["read_id"] = df["qseqid"]
        rows.append(df)
    if not rows:
        return pd.DataFrame(columns=["read_id", "family", "bitscore", "evalue", "sample"])
    return pd.concat(rows, ignore_index=True)


# orfm's ORF headers are exactly "<original read ID>_<orf>_<frame>_<n>"
# (verified against real orfm output: e.g. ".../1_1_1_1", ".../1_2_5_5"),
# so stripping the trailing three underscore-separated numbers after the
# mate marker recovers the original read ID it was called from.
ORF_SUFFIX_RE = re.compile(r"^(.*/[12])_\d+_\d+_\d+$")


def load_hmm_calls(hmm_dir: Path) -> pd.DataFrame:
    """Parse hmmscan --tblout files into a long dataframe of read -> family calls."""
    rows = []
    for tblout in hmm_dir.glob("*.hmmscan.tblout"):
        sample = tblout.stem.replace(".hmmscan", "")
        with open(tblout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                target_hmm_name = fields[0]   # family name (HMM model name)
                query_name = fields[2]        # translated ORF id (orfm output)
                bitscore = float(fields[5])
                m = ORF_SUFFIX_RE.match(query_name)
                read_id = m.group(1) if m else query_name
                rows.append({
                    "sample": sample,
                    "read_id": read_id,
                    "family": target_hmm_name,
                    "bitscore": bitscore,
                })
    return pd.DataFrame(rows)


def best_hit_per_read(calls: pd.DataFrame) -> pd.DataFrame:
    """
    Reduce a long table of (read, family, score) hits down to a single
    best-scoring family call per read (per sample). Without this, a read
    that weakly matches several unrelated HMMs/DIAMOND targets gets counted
    as a false positive once per spurious hit instead of once per read,
    which badly inflates the false-positive count for promiscuous families
    (e.g. ectB / betL, whose Pfam domains are shared across many unrelated
    protein families).
    """
    if calls.empty:
        return calls
    idx = calls.groupby(["sample", "read_id"])["bitscore"].idxmax()
    return calls.loc[idx].reset_index(drop=True)


def label_reads(calls: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """
    Join calls to per-read truth on read_id. Uses a left join FROM truth so
    that reads with zero hits at all (never appearing in `calls`) still show
    up with family_called="NO_CALL" and get scored as misses, instead of
    silently vanishing from the metric the way an inner join would.

    A read only counts as a genuine positive if it's both labeled "positive"
    AND actually overlaps the true CDS (see read_truth.tsv/overlaps_cds) --
    a read from a positive construct that lands entirely in flanking DNA
    carries no gene signal and must not be scored as something the caller
    was supposed to recall.
    """
    calls = best_hit_per_read(calls)
    merged = truth.merge(calls, on="read_id", how="left", suffixes=("_true", "_called"))
    merged["family_called"] = merged["family_called"].fillna("NO_CALL")
    # sklearn's roc_curve/precision_recall_curve reject +/-inf, so use a
    # sentinel below any real bitscore instead (uncalled reads must never
    # pass a positive score threshold)
    min_real_score = merged["bitscore"].min()
    sentinel = (min_real_score - 1.0) if pd.notna(min_real_score) else -1.0
    merged["bitscore"] = merged["bitscore"].fillna(sentinel)
    merged["is_genuine_positive"] = (
        (merged["label"] == "positive") & merged["overlaps_cds"].astype(bool)
    )
    merged["correct"] = merged["is_genuine_positive"] & (
        merged["family_called"] == merged["family_true"]
    )
    return merged


def plot_curves(merged: pd.DataFrame, method: str, score_col: str, out_dir: Path,
                 higher_is_better: bool = True) -> dict:
    y_true = merged["correct"].astype(int)
    scores = merged[score_col] if higher_is_better else -merged[score_col]

    if y_true.nunique() < 2:
        return {}

    fpr, tpr, _ = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(y_true, scores)
    pr_auc = auc(rec, prec)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(fpr, tpr, label=f"AUC={roc_auc:.3f}")
    axes[0].plot([0, 1], [0, 1], "k--", lw=0.5)
    axes[0].set_xlabel("False positive rate")
    axes[0].set_ylabel("True positive rate")
    axes[0].set_title(f"{method} ROC")
    axes[0].legend()

    axes[1].plot(rec, prec, label=f"AUC={pr_auc:.3f}")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title(f"{method} PR")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(out_dir / f"{method}_roc_pr.png", dpi=150)
    plt.close(fig)

    return {"roc_auc": roc_auc, "pr_auc": pr_auc}


def summarize_best_threshold(merged: pd.DataFrame, method: str, score_col: str,
                              families: list[str]) -> list[dict]:
    """
    For each family, sweep the score column and report the threshold that
    maximizes F1 for calls made to that family. This shows the best
    achievable precision/recall if a per-family score cutoff were applied
    on top of the best-hit-per-read call (e.g. a stricter HMM bit-score
    cutoff, or a stricter DIAMOND bitscore/e-value cutoff), instead of just
    accepting every best-hit call regardless of score.
    """
    rows = []
    for family in families:
        fam_calls = merged[merged["family_called"] == family]
        if fam_calls.empty:
            continue
        is_true_positive_source = (fam_calls["family_true"] == family) & fam_calls["is_genuine_positive"]
        scores = fam_calls[score_col].to_numpy()
        y_true = is_true_positive_source.to_numpy().astype(int)
        if y_true.sum() == 0 or y_true.sum() == len(y_true):
            continue

        best = {"f1": -1.0}
        for thr in sorted(set(scores)):
            keep = scores >= thr
            tp = int((y_true[keep] == 1).sum())
            fp = int((y_true[keep] == 0).sum())
            fn = int((y_true == 1).sum() - tp)
            precision = tp / (tp + fp) if (tp + fp) else float("nan")
            recall = tp / (tp + fn) if (tp + fn) else float("nan")
            f1 = (2 * precision * recall / (precision + recall)
                  if precision == precision and recall == recall and (precision + recall) else 0.0)
            if f1 > best["f1"]:
                best = {"f1": f1, "threshold": thr, "tp": tp, "fp": fp, "fn": fn,
                        "precision": precision, "recall": recall}

        rows.append({
            "method": method, "family": family,
            "best_threshold": round(best["threshold"], 2),
            "tp": best["tp"], "fp": best["fp"], "fn": best["fn"],
            "precision": round(best["precision"], 3) if best["precision"] == best["precision"] else "",
            "recall": round(best["recall"], 3) if best["recall"] == best["recall"] else "",
            "f1": round(best["f1"], 3),
        })
    return rows


def summarize(merged: pd.DataFrame, method: str, families: list[str]) -> list[dict]:
    rows = []
    for family in families:
        fam_calls = merged[merged["family_true"] == family]
        if fam_calls.empty:
            continue
        tp = ((fam_calls["family_called"] == family) & fam_calls["is_genuine_positive"]).sum()
        fn = (fam_calls["is_genuine_positive"] & (fam_calls["family_called"] != family)).sum()
        fp = ((fam_calls["family_called"] == family) & ~fam_calls["is_genuine_positive"]).sum()
        precision = tp / (tp + fp) if (tp + fp) else float("nan")
        recall = tp / (tp + fn) if (tp + fn) else float("nan")
        f1 = (2 * precision * recall / (precision + recall)
              if precision and recall and (precision + recall) else float("nan"))
        rows.append({
            "method": method, "family": family, "tp": tp, "fp": fp, "fn": fn,
            "precision": round(precision, 3) if precision == precision else "",
            "recall": round(recall, 3) if recall == recall else "",
            "f1": round(f1, 3) if f1 == f1 else "",
        })
    return rows


def summarize_by_length(merged: pd.DataFrame, method: str, families: list[str],
                         bucket_size: int = 50) -> list[dict]:
    """Same precision/recall/F1 as summarize(), but broken out per read-length
    bucket, so a method that only works on long reads doesn't get credit it
    wouldn't earn on a typical short-read metagenomic sample."""
    rows = []
    buckets = (merged["read_length"] // bucket_size * bucket_size).astype(int)
    for bucket in sorted(buckets.unique()):
        bucket_rows = merged[buckets == bucket]
        for row in summarize(bucket_rows, method, families):
            row["read_length_bucket"] = f"{bucket}-{bucket + bucket_size - 1}bp"
            rows.append(row)
    return rows


def summarize_background_fpr(merged: pd.DataFrame, method: str) -> dict | None:
    """Real background reads (family == 'background') carry no osmoadaptation
    gene content at all -- any call made on one is a false positive, and
    this can't be folded into summarize()'s per-family tables since
    'background' isn't a real gene family."""
    bg = merged[merged["family_true"] == "background"]
    if bg.empty:
        return None
    n_reads = len(bg)
    n_false_positive = (bg["family_called"] != "NO_CALL").sum()
    return {
        "method": method,
        "n_background_reads": n_reads,
        "n_false_positive_calls": int(n_false_positive),
        "background_fpr": round(n_false_positive / n_reads, 4),
    }


def build_profile_cascade_config(
    summary_rows: list[dict], best_threshold_rows: list[dict],
    precision_threshold: float,
) -> list[dict]:
    """
    For each family, pair DIAMOND's benchmark precision (accept-every-
    best-hit, the same operating point osmotool's `profile` actually uses
    in production) with HMM's best short-read bitscore threshold from
    summarize_best_threshold(). Families below precision_threshold are
    flagged for osmotool's profile-mode DIAMOND+HMM cascade: DIAMOND stays
    the primary caller for speed, but reads it assigns to a flagged family
    get a second opinion from hmmscan (raw bitscore, NOT --cut_ga -- these
    are short reads, not full-length ORFs, so the cascade needs HMM's
    short-read-calibrated threshold here, not the GA cutoff used in
    `annotate`) before the call is kept.

    Families with perfect/near-perfect DIAMOND precision are left
    unflagged so the cascade only costs extra hmmscan time where DIAMOND
    has a demonstrated real specificity problem, not on every read.
    """
    diamond_precision = {
        r["family"]: r["precision"] for r in summary_rows if r["method"] == "diamond"
    }
    hmm_threshold = {
        r["family"]: r["best_threshold"] for r in best_threshold_rows if r["method"] == "hmm"
    }

    rows = []
    for family, precision in diamond_precision.items():
        precision_val = precision if precision != "" else float("nan")
        needs_cascade = (
            precision_val == precision_val  # not NaN
            and precision_val < precision_threshold
            and family in hmm_threshold
        )
        rows.append({
            "family": family,
            "diamond_precision": precision_val,
            "needs_cascade_check": "yes" if needs_cascade else "no",
            "hmm_short_read_threshold": hmm_threshold.get(family, "") if needs_cascade else "",
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--cascade-out", type=Path, default=None,
                     help="Where to write profile_cascade.tsv (osmotool's "
                          "profile-mode DIAMOND+HMM cascade config). "
                          "Defaults to <out>/profile_cascade.tsv.")
    ap.add_argument("--cascade-precision-threshold", type=float, default=0.95,
                     help="Families with DIAMOND precision below this are "
                          "flagged for the profile-mode cascade check.")
    args = ap.parse_args()

    families = load_family_names(args.families)
    args.out.mkdir(parents=True, exist_ok=True)
    truth = load_read_truth(args.results / "reads")

    diamond_calls = load_diamond_calls(args.results / "diamond")
    hmm_calls = load_hmm_calls(args.results / "hmm")

    summary_rows = []
    best_threshold_rows = []
    length_stratified_rows = []
    background_fpr_rows = []
    all_curve_stats = {}

    if not diamond_calls.empty:
        diamond_merged = label_reads(diamond_calls, truth)
        summary_rows += summarize(diamond_merged, "diamond", families)
        best_threshold_rows += summarize_best_threshold(diamond_merged, "diamond", "bitscore", families)
        length_stratified_rows += summarize_by_length(diamond_merged, "diamond", families)
        bg_row = summarize_background_fpr(diamond_merged, "diamond")
        if bg_row:
            background_fpr_rows.append(bg_row)
        all_curve_stats["diamond"] = plot_curves(
            diamond_merged, "diamond", "bitscore", args.out, higher_is_better=True)
    else:
        print("WARNING: no DIAMOND alignment files found in results/diamond")

    if not hmm_calls.empty:
        hmm_merged = label_reads(hmm_calls, truth)
        summary_rows += summarize(hmm_merged, "hmm", families)
        best_threshold_rows += summarize_best_threshold(hmm_merged, "hmm", "bitscore", families)
        length_stratified_rows += summarize_by_length(hmm_merged, "hmm", families)
        bg_row = summarize_background_fpr(hmm_merged, "hmm")
        if bg_row:
            background_fpr_rows.append(bg_row)
        all_curve_stats["hmm"] = plot_curves(
            hmm_merged, "hmm", "bitscore", args.out, higher_is_better=True)
    else:
        print("WARNING: no HMM tblout files found in results/hmm")

    summary_path = args.out / "summary.tsv"
    if summary_rows:
        with open(summary_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"Per-family precision/recall/F1 written to {summary_path}")

    best_threshold_path = args.out / "best_threshold_summary.tsv"
    if best_threshold_rows:
        with open(best_threshold_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(best_threshold_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(best_threshold_rows)
        print(f"Per-family best-F1 score-threshold summary written to {best_threshold_path}")

    length_path = args.out / "summary_by_read_length.tsv"
    if length_stratified_rows:
        with open(length_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(length_stratified_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(length_stratified_rows)
        print(f"Read-length-stratified precision/recall/F1 written to {length_path}")

    background_path = args.out / "background_fpr.tsv"
    if background_fpr_rows:
        with open(background_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(background_fpr_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(background_fpr_rows)
        print(f"Background false-positive rate written to {background_path}")

    if all_curve_stats:
        print("\nOverall ROC/PR AUC by method:")
        for method, stats in all_curve_stats.items():
            print(f"  {method}: {stats}")

    cascade_rows = build_profile_cascade_config(
        summary_rows, best_threshold_rows, args.cascade_precision_threshold)
    if cascade_rows:
        cascade_path = args.cascade_out or (args.out / "profile_cascade.tsv")
        with open(cascade_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(cascade_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(cascade_rows)
        n_flagged = sum(1 for r in cascade_rows if r["needs_cascade_check"] == "yes")
        print(f"\nProfile-mode DIAMOND+HMM cascade config written to {cascade_path} "
              f"({n_flagged}/{len(cascade_rows)} families flagged, "
              f"precision threshold {args.cascade_precision_threshold})")


if __name__ == "__main__":
    main()
