#!/usr/bin/env python3
"""
11_compute_metrics.py — score DIAMOND vs. HMM against the simulated-read truth table

Parses:
  - results/diamond/<sample>.blastx.tsv (osmotool --keep_aln output)
  - results/hmm/<sample>.hmmscan.tblout (hmmscan output, raw bit scores)
  - results/reads/truth.tsv (from 09_simulate_reads.py)

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


def load_truth(reads_dir: Path) -> pd.DataFrame:
    truth = pd.read_csv(reads_dir / "truth.tsv", sep="\t")
    # read_id prefix in simulated FASTQ typically embeds the source header;
    # downstream parsing keys on the UniProt ID portion of source_header.
    truth["uniprot_id"] = truth["source_header"].str.split("|").str[1]
    return truth


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
        # src/osmotool/functions.py (NOT the 12-column BLAST default)
        df.columns = [
            "qseqid", "sseqid", "pident", "qcovhsp", "length", "evalue", "bitscore",
        ][: df.shape[1]]
        # keep only the best (highest-bitscore) hit per read for scoring
        df = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
        df["sample"] = sample
        df["family"] = df["sseqid"].str.split("|").str[0]
        rows.append(df)
    if not rows:
        return pd.DataFrame(columns=["qseqid", "family", "bitscore", "evalue", "sample"])
    return pd.concat(rows, ignore_index=True)


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
                query_name = fields[2]        # translated ORF / read id
                bitscore = float(fields[5])
                rows.append({
                    "sample": sample,
                    "qseqid": query_name,
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
    idx = calls.groupby(["sample", "qseqid"])["bitscore"].idxmax()
    return calls.loc[idx].reset_index(drop=True)


def label_reads(calls: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """
    Join calls to truth via UniProt ID embedded in the read name (simulators
    are expected to preserve the source contig/sequence header in read IDs;
    adjust the extraction regex here if your simulator's naming differs).
    Only the single best-scoring hit per read is kept (see
    best_hit_per_read) so that one read can only ever be one false positive,
    not one-per-spurious-model-hit.
    """
    calls = best_hit_per_read(calls)
    calls = calls.copy()
    calls["uniprot_id"] = calls["qseqid"].str.extract(r"\|([A-Z0-9]+)\|")[0]
    merged = calls.merge(truth[["uniprot_id", "family", "label"]],
                          on="uniprot_id", suffixes=("_called", "_true"))
    merged["correct"] = (
        (merged["family_called"] == merged["family_true"])
        & (merged["label"] == "positive")
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
        is_true_positive_source = (fam_calls["family_true"] == family) & (fam_calls["label"] == "positive")
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
        tp = ((fam_calls["family_called"] == family) & (fam_calls["label"] == "positive")).sum()
        fn = ((fam_calls["label"] == "positive") & (fam_calls["family_called"] != family)).sum()
        fp = ((fam_calls["family_called"] == family) & (fam_calls["label"] != "positive")).sum()
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    args = ap.parse_args()

    families = load_family_names(args.families)
    args.out.mkdir(parents=True, exist_ok=True)
    truth = load_truth(args.results / "reads")

    diamond_calls = load_diamond_calls(args.results / "diamond")
    hmm_calls = load_hmm_calls(args.results / "hmm")

    summary_rows = []
    best_threshold_rows = []
    all_curve_stats = {}

    if not diamond_calls.empty:
        diamond_merged = label_reads(diamond_calls, truth)
        summary_rows += summarize(diamond_merged, "diamond", families)
        best_threshold_rows += summarize_best_threshold(diamond_merged, "diamond", "bitscore", families)
        all_curve_stats["diamond"] = plot_curves(
            diamond_merged, "diamond", "bitscore", args.out, higher_is_better=True)
    else:
        print("WARNING: no DIAMOND alignment files found in results/diamond")

    if not hmm_calls.empty:
        hmm_merged = label_reads(hmm_calls, truth)
        summary_rows += summarize(hmm_merged, "hmm", families)
        best_threshold_rows += summarize_best_threshold(hmm_merged, "hmm", "bitscore", families)
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

    if all_curve_stats:
        print("\nOverall ROC/PR AUC by method:")
        for method, stats in all_curve_stats.items():
            print(f"  {method}: {stats}")


if __name__ == "__main__":
    main()
