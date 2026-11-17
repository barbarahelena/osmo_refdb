#!/usr/bin/env python3
"""
compute_metrics.py — score DIAMOND vs. HMM against the simulated-read truth table

Parses:
  - results/diamond/<sample>.blastx.tsv (osmotool --keep_aln output)
  - results/hmm/<sample>.hmmscan.tblout (hmmscan output, raw bit scores)
  - results/reads/read_truth.tsv (from compute_read_truth.py; per-read
    truth based on actual CDS overlap, not just source-construct label)

For each family and method, computes read-level precision/recall/F1. A read
is reduced to its single best-scoring family call (see best_hit_per_read) so
that one read can only ever count as one false positive, not once per
spurious model hit -- important for promiscuous domains (e.g. ectB, betL)
that legitimately score weakly against many unrelated HMMs/DIAMOND targets.

False positives are counted three ways for family F (summary.tsv):
  fp            reads from F's OWN constructs called F: its held-out hard
                negatives (same Pfam domain, different gene) and flanking
                DNA. precision / f1 use only these -- the definition every
                release up to v8 reported, kept for comparison.
  fp_cross      genuine gene reads of any OTHER family called F: real
                confusion between panel genes (e.g. a ktrB read called
                trkH). Before 2026-11-04 these were silently dropped: a
                misassigned read only counted as a miss for its own family.
                precision_panel / f1_panel count fp + fp_cross. A
                diagnostic, not a decision input: close paralogs in the
                panel (opu ATPases/permeases, trkH/ktrB/ktrD, ...) are
                expected to confuse each other, and whether that matters
                depends on whether a result is read per gene or per system.
  fp_other_negatives
                OTHER families' hard-negative and flank reads called F.
                Reported, but in no precision: sibling families share one
                negative_query (trkH/ktrB/ktrD, proX/opuAC/opuBC/opuCC,
                trkA/ktrA, gshA/gshF), so their hard negatives are further
                draws from F's own hard-negative pool, and counting them
                would weight F's hard negatives by how many siblings it
                has. On v8 they were ~99% of all cross-family calls (trkH:
                12,097 of 12,191), real-gene confusion ~1%.
  fp_background reads from the shared random bacterial background
                (build_random_background.py, family "random_bg") called
                F, also as background_fp_per_million. Kept out of every
                precision: how many unrelated reads to mix in is arbitrary,
                so the rate per million unrelated reads is the number to
                read (it scales directly to a real metagenome).

  fp_all        fp + fp_cross + fp_other_negatives: every non-background
                read called F, whatever construct it came from.
                precision_all / f1_all use it. This is the denominator to
                compare gene level with group level: a group's fp counts
                every member's hard negatives and flanks, so its precision
                charges a sibling's negatives that per-gene precision never
                sees. With fp_all both levels count the same reads: summed
                over a group's members, fp_all is the group's fp_all plus
                its within_group_swaps, which the group instead counts as
                tp. The report-level and cascade decision (below) uses f1_all.

summary_by_group.tsv repeats these metrics per families.yaml
`paralog_group`, scoring a group's members as one unit: a read called as
any member of its true family's group is a hit (within_group_swaps counts
how many of those were gene-level misses).

Calibration mode (no --cascade-config; run on the calibration reads with
run_benchmark.sh CASCADE=score) decides, per reporting unit, the report
level and the cascade together, and writes both release manifests:
osmo_refdb.profile_cascade.tsv (which families get the HMM second opinion)
and osmo_refdb.profile_report_levels.tsv (what osmotool profile reports, per
unit, with the calibration numbers each choice rests on). See decide_units.
Benchmark mode (--cascade-config, --report-levels; the test reads) only
scores those choices, adding summary_by_report_unit.tsv.

Precision depends on how many negative reads are mixed in, which is set by
fetch caps and split fractions, not biology. hard_negative_fp_per_million
(hard-negative reads called as the family, per million hard-negative reads)
and background_fp_per_million don't: read those two side by side to compare
families, groups or releases. At group level the hard-negative rate is the
chance a hard-negative read lands on ANY member, so it also catches negatives
"claimed" by a sibling, which per-gene precision never counts (a negative
from opuBA's sample called opuAA is neither opuBA's error nor opuAA's own
negative).

Also sweeps DIAMOND's/HMM's bit-score thresholds to produce ROC/PR-style
curves (saved as PNGs) plus two summary tables:
  - summary.tsv: precision/recall/F1 accepting every best-hit call as-is
  - best_threshold_summary.tsv: precision/recall/F1 at the per-family score
    threshold that maximizes F1 (i.e. the best achievable operating point if
    a stricter score cutoff were applied on top of the best-hit call)

Usage:
  python compute_metrics.py --results results/ --out results/metrics --families families.yaml
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


RANDOM_BG_FAMILY = "random_bg"             # must match 01f / 09
REAL_BG_FAMILY = "background"              # --background real metagenome (09)
BACKGROUND_SOURCES = (RANDOM_BG_FAMILY, REAL_BG_FAMILY)


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def load_read_truth(reads_dir: Path) -> pd.DataFrame:
    """
    Load per-read ground truth from compute_read_truth.py's read_truth.tsv
    (read_id, family, label, uniprot_id, overlaps_cds, overlap_bp). Unlike the
    contig-level truth.tsv (one row per source protein), this reflects whether
    each individual simulated read actually overlaps the true CDS -- a read
    from a "positive" construct that lands entirely in the random flanking
    DNA carries no gene signal and must not be scored as a true positive.
    """
    return pd.read_csv(reads_dir / "read_truth.tsv", sep="\t")


MATE_TAG_RE = re.compile(r"^(.*/[12])_R[12]$")


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
        # same tie-break as osmotool's quantifier.select_best_hits: higher
        # bitscore, then lower evalue, then lexicographically smallest sseqid
        df = (df.sort_values(["bitscore", "evalue", "sseqid"],
                             ascending=[False, True, True], kind="mergesort")
                .drop_duplicates("qseqid"))
        df["sample"] = sample
        df["family"] = df["sseqid"].str.split("|").str[0]
        # osmotool >= ab3fb9e appends _R1/_R2 to paired read IDs (so two
        # mates can't collide); the simulated reads' own IDs already end in
        # /1 or /2, so stripping that tag recovers the read_truth.tsv ID.
        # Older osmotool reports the raw read ID, which this leaves unchanged.
        df["read_id"] = df["qseqid"].str.replace(MATE_TAG_RE, r"\1", regex=True)
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


CASCADE_ORF_SUFFIX_RE = re.compile(r"^(.*)_\d+_\d+_\d+$")   # osmotool.functions.ORF_SUFFIX_RE


def load_cascade_config(path: Path) -> dict[str, float]:
    """{family: hmm_short_read_threshold} for needs_cascade_check == yes;
    same reading as osmotool.functions.load_cascade_config."""
    cutoffs: dict[str, float] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["needs_cascade_check"] == "yes":
                cutoffs[row["family"]] = float(row["hmm_short_read_threshold"])
    return cutoffs


def apply_profile_cascade(diamond_calls: pd.DataFrame, diamond_dir: Path,
                          cutoffs: dict[str, float]) -> pd.DataFrame:
    """
    Per-read calls of `osmotool profile` with its default DIAMOND+HMM cascade.

    osmotool reports only per-family counts, so the per-read outcome is
    rebuilt from what --keep_aln leaves behind: the raw DIAMOND alignments
    (diamond_calls, one best hit per read) and osmotool's own
    <sample>.cascade.hmmsearch.tblout. The rule is
    quantifier.apply_cascade_check's: a read called to a flagged family is
    kept only if its best raw hmmsearch bitscore against THAT family's HMM
    reaches the family's threshold; reads on unflagged families pass
    untouched. verify_cascade_counts() checks the result against the
    gene_counts.tsv osmotool wrote in the same run.
    """
    if diamond_calls.empty or not cutoffs:
        return diamond_calls.copy()
    score_rows = []
    for tblout in diamond_dir.glob("*.cascade.hmmsearch.tblout"):
        sample = tblout.name.replace(".cascade.hmmsearch.tblout", "")
        with open(tblout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                m = CASCADE_ORF_SUFFIX_RE.match(fields[0])
                score_rows.append((sample, m.group(1) if m else fields[0],
                                   fields[2], float(fields[5])))
    scores = pd.DataFrame(score_rows, columns=["sample", "read_id", "family", "hmm_score"])
    scores = scores.groupby(["sample", "read_id", "family"], as_index=False)["hmm_score"].max()
    calls = diamond_calls.merge(scores, on=["sample", "read_id", "family"], how="left")
    calls["hmm_score"] = calls["hmm_score"].fillna(0.0)
    threshold = calls["family"].map(cutoffs)
    keep = threshold.isna() | (calls["hmm_score"] >= threshold)
    return calls[keep].drop(columns="hmm_score").reset_index(drop=True)


def verify_cascade_counts(cascade_calls: pd.DataFrame, diamond_dir: Path) -> None:
    """Compare the rebuilt per-read cascade calls with osmotool's own
    gene_counts.tsv (each paired read counts 0.5). Families osmotool drops
    from its report (profile_excluded_families, decoys) are simply absent
    there and are skipped. Prints the number of disagreeing cells."""
    rebuilt = (cascade_calls.groupby(["sample", "family"]).size() * 0.5).to_dict()
    n_cells = n_bad = 0
    worst: list[tuple] = []
    for counts_file in diamond_dir.glob("*.gene_counts.tsv"):
        sample = counts_file.name.replace(".gene_counts.tsv", "")
        df = pd.read_csv(counts_file, sep="\t", comment="#")
        for family, raw in zip(df["gene"], df["raw_count"]):
            n_cells += 1
            mine = rebuilt.get((sample, family), 0.0)
            if abs(mine - raw) > 1e-6:
                n_bad += 1
                worst.append((abs(mine - raw), sample, family, raw, mine))
    print(f"Cascade check: rebuilt per-read calls vs osmotool gene_counts.tsv: "
          f"{n_cells - n_bad}/{n_cells} sample x family cells agree")
    for _, sample, family, raw, mine in sorted(worst, reverse=True)[:10]:
        print(f"  MISMATCH {sample} {family}: osmotool {raw}, rebuilt {mine}")


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
    if not calls.empty:
        # A read whose best hit is a decoy (DIAMOND "<family>_decoy" refs,
        # 08a; HMM decoy models, 05b) was won away from every real family:
        # osmotool drops these labels from output, so score it as no call.
        # Done after best_hit_per_read so the decoy still competes.
        calls = calls[~calls["family"].astype(str).str.endswith("_decoy")]
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
    # Masks _metrics() needs for every family/group/length bucket, computed
    # once here instead of rescanning the string columns each call.
    merged["_bg_source"] = merged["family_true"].isin(BACKGROUND_SOURCES)
    merged["_random_bg"] = merged["family_true"] == RANDOM_BG_FAMILY
    merged["_hard_neg"] = merged["label"] == "negative"
    merged["family_true"] = merged["family_true"].astype("category")
    merged["family_called"] = merged["family_called"].astype("category")
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
    # Background reads are left out of the sweep: their volume is arbitrary,
    # and including them would let it move the threshold the profile-mode
    # cascade ships with. Cross-family reads stay in (they were always
    # counted here, unlike in summary.tsv before fp_cross existed).
    merged = merged[~merged["family_true"].isin(BACKGROUND_SOURCES)]
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


def _ratio(num: int, den: int) -> float:
    return num / den if den else float("nan")


def _f1(precision: float, recall: float) -> float:
    if precision != precision or recall != recall or not (precision + recall):
        return float("nan")
    return 2 * precision * recall / (precision + recall)


def _r3(x: float):
    return round(x, 3) if x == x else ""


def _metrics(merged: pd.DataFrame, members: set[str], n_random_bg: int) -> dict | None:
    """Counts and metrics for a set of families scored as one unit: a single
    family (summary.tsv) or a paralog group (summary_by_group.tsv). A read
    from a member's construct called as any member counts as a hit, so at
    group level confusion between siblings is not an error. See the module
    docstring for the false-positive columns."""
    own = merged["family_true"].isin(members)
    if not own.any():
        return None
    called = merged["family_called"].isin(members)
    genuine = merged["is_genuine_positive"]
    other = ~own & ~merged["_bg_source"] & called
    tp = int((own & called & genuine).sum())
    fn = int((own & genuine & ~called).sum())
    fp = int((own & called & ~genuine).sum())
    fp_cross = int((other & genuine).sum())
    fp_other_neg = int((other & ~genuine).sum())
    fp_bg = int((merged["_random_bg"] & called).sum())
    hard_neg = own & merged["_hard_neg"]
    n_hard_neg = int(hard_neg.sum())
    fp_hard_neg = int((hard_neg & called).sum())
    precision = _ratio(tp, tp + fp)
    recall = _ratio(tp, tp + fn)
    precision_panel = _ratio(tp, tp + fp + fp_cross)
    fp_all = fp + fp_cross + fp_other_neg
    precision_all = _ratio(tp, tp + fp_all)
    return {
        "tp": tp, "fp": fp, "fn": fn,
        "precision": _r3(precision),
        "recall": _r3(recall),
        "f1": _r3(_f1(precision, recall)),
        "fp_cross": fp_cross,
        "fp_other_negatives": fp_other_neg,
        "precision_panel": _r3(precision_panel),
        "f1_panel": _r3(_f1(precision_panel, recall)),
        "fp_all": fp_all,
        "precision_all": _r3(precision_all),
        "f1_all": _r3(_f1(precision_all, recall)),
        "n_hard_negative_reads": n_hard_neg,
        "hard_negative_fp_per_million": (round(fp_hard_neg / n_hard_neg * 1e6, 1)
                                         if n_hard_neg else ""),
        "fp_background": fp_bg if n_random_bg else "",
        "background_fp_per_million": (round(fp_bg / n_random_bg * 1e6, 1)
                                      if n_random_bg else ""),
        # genuine member reads called as a DIFFERENT member: a miss at gene
        # level, a hit at group level (always 0 for a single family)
        "within_group_swaps": (int((own & genuine & called & (
            merged["family_called"].astype(str) != merged["family_true"].astype(str))).sum())
            if len(members) > 1 else 0),
    }


def summarize(merged: pd.DataFrame, method: str, families: list[str],
              groups: dict[str, str] | None = None) -> list[dict]:
    """Per-family counts and metrics; see the module docstring for the three
    false-positive columns."""
    groups = groups or {}
    n_random_bg = int(merged["_random_bg"].sum())
    rows = []
    for family in families:
        m = _metrics(merged, {family}, n_random_bg)
        if m is None:
            continue
        del m["within_group_swaps"]
        rows.append({"method": method, "family": family, **m,
                     "paralog_group": groups.get(family, "")})
    return rows


def summarize_groups(merged: pd.DataFrame, method: str, families: list[str],
                     groups: dict[str, str]) -> list[dict]:
    """Per paralog group (families.yaml paralog_group): the same metrics as
    summary.tsv with the group's members scored as one unit. Next to each
    row, the members' summed per-gene tp so the gain from not counting
    sibling swaps is visible."""
    n_random_bg = int(merged["_random_bg"].sum())
    members_by_group: dict[str, list[str]] = {}
    for family in families:
        if family in groups:
            members_by_group.setdefault(groups[family], []).append(family)
    rows = []
    for group, members in members_by_group.items():
        m = _metrics(merged, set(members), n_random_bg)
        if m is None:
            continue
        rows.append({"method": method, "paralog_group": group,
                     "members": ",".join(members), **m})
    return rows


REPORT_LEVELS = ("gene", "group", "none")


def load_family_meta(path: Path) -> dict[str, dict]:
    """{family: {"group", "report_level", "scope"}} from families.yaml, in
    file order. report_level is checked here: one of gene/group/none,
    "group" only for a family with a paralog_group, and "none" for
    scope: annotate_only (those are
    excluded from profile output whatever calibration says). A group's
    members either all say "group" (reported as one unit) or each say gene
    or none on their own. An unset
    report_level means "take the calibration rule's level" (markers.yaml,
    subset files)."""
    with open(path) as fh:
        data = yaml.safe_load(fh)
    meta: dict[str, dict] = {}
    for fam in data["families"]:
        meta[fam["name"]] = {"group": fam.get("paralog_group") or "",
                             "report_level": fam.get("report_level") or "",
                             "scope": fam.get("scope") or ""}
    errors = []
    levels_by_group: dict[str, set[str]] = {}
    for name, m in meta.items():
        level = m["report_level"]
        if level and level not in REPORT_LEVELS:
            errors.append(f"{name}: report_level '{level}' is not one of {', '.join(REPORT_LEVELS)}")
        if level == "group" and not m["group"]:
            errors.append(f"{name}: report_level group needs a paralog_group")
        if m["scope"] == "annotate_only" and level not in ("", "none"):
            errors.append(f"{name}: scope annotate_only needs report_level none, not {level}")
        if m["group"]:
            levels_by_group.setdefault(m["group"], set()).add(level)
    for group, levels in levels_by_group.items():
        if "group" in levels and len(levels) > 1:
            errors.append(f"paralog_group {group}: report_level group must be set on every "
                          f"member or none ({', '.join(sorted(lv or 'unset' for lv in levels))})")
    if errors:
        raise SystemExit("families.yaml report_level:\n  " + "\n  ".join(errors))
    return meta


def _f1_value(row: dict | None) -> float:
    """f1_all of a summarize()/summarize_groups() row as a number; a unit
    with no calibration reads or no calls at all counts as 0."""
    if row is None or row["f1_all"] == "":
        return 0.0
    return float(row["f1_all"])


def _best_cascade(off: dict | None, on: dict | None) -> tuple[bool, float]:
    """(cascade on?, f1_all) for the better of the two options. The cascade
    costs hmmsearch time and recall, so it is chosen only if it raises
    f1_all (on the 3-decimal values in the manifest, so anyone can redo the
    comparison from it)."""
    f_off = _f1_value(off)
    if on is None:
        return False, f_off
    f_on = _f1_value(on)
    return (True, f_on) if f_on > f_off else (False, f_off)


def _fmt(x) -> str:
    return "NA" if x is None or x == "" else str(x)


def decide_units(meta: dict[str, dict], thresholds: dict[str, float],
                 gene_off: dict[str, dict], gene_on: dict[str, dict],
                 group_off: dict[str, dict], group_on: dict[str, dict],
                 min_f1: float) -> tuple[list[dict], dict[str, dict]]:
    """
    Decide report level and cascade together, per family or paralog group
    (families.yaml paralog_group), on the calibration reads.

    gene_off / gene_on: {family: summarize() row} on the calibration reads
    without the cascade and with every family that has an HMM threshold
    cascaded. group_off / group_on: {group: summarize_groups() row}, same.
    The cascade drops reads but never moves one to another family, so a
    family's numbers depend only on its own flag and a group's only on its
    members': these four tables score every option.

    Rule, on f1_all (see the module docstring), each level taking the
    better of its cascade options:
      a family without a group: gene if it reaches min_f1, else none;
      a paralog group: gene for every member if every member reaches
        min_f1; else group (one unit, one cascade option for all members)
        if the group reaches min_f1; else each member on its own, gene if
        it reaches min_f1, none if not.
    families.yaml's report_level is what gets reported; the rule's level is
    written next to it (calibration_level) and a mismatch is printed as a
    warning. The cascade follows the reported level (off for none).

    Returns (manifest rows: one per reporting unit, i.e. a paralog group on
    "group" and a family otherwise; {family: cascade-config row}).
    """
    blocks: list[tuple[list[str], str]] = []
    seen_groups: set[str] = set()
    for family, m in meta.items():
        if not m["group"]:
            blocks.append(([family], ""))
        elif m["group"] not in seen_groups:
            seen_groups.add(m["group"])
            blocks.append(([f for f, mm in meta.items() if mm["group"] == m["group"]], m["group"]))

    def f1_cell(table: dict, key: str) -> str:
        return _fmt(table.get(key, {}).get("f1_all"))

    def chosen(row: dict | None) -> dict:
        row = row or {}
        return {k: row.get(k, "") for k in ("tp", "fp_all", "fn", "precision_all", "recall", "f1_all")}

    manifest: list[dict] = []
    cascade: dict[str, dict] = {}
    for members, group in blocks:
        gene_choice = {f: _best_cascade(gene_off.get(f), gene_on.get(f) if f in thresholds else None)
                       for f in members}
        group_has_on = any(f in thresholds for f in members)
        group_choice = (_best_cascade(group_off.get(group), group_on.get(group) if group_has_on else None)
                        if group else None)

        weak = [f for f in members if gene_choice[f][1] < min_f1]
        calib: dict[str, tuple[str, str]] = {}
        for f in members:
            f1 = gene_choice[f][1]
            if not weak:
                calib[f] = ("gene", f"every member of {group} reaches f1_all {min_f1} at gene level ({f} {f1:.3f})"
                            if group else f"f1_all {f1:.3f} >= {min_f1}")
            elif group_choice[1] >= min_f1 if group else False:
                calib[f] = ("group", f"{', '.join(weak)} below f1_all {min_f1} at gene level; "
                                     f"group f1_all {group_choice[1]:.3f} >= {min_f1}")
            else:
                tail = f"; group {group} f1_all {group_choice[1]:.3f} < {min_f1}" if group else ""
                calib[f] = (("gene", f"f1_all {f1:.3f} >= {min_f1} at gene level{tail}") if f not in weak
                            else ("none", f"f1_all {f1:.3f} < {min_f1} at gene level{tail}"))

        level: dict[str, str] = {}
        source: dict[str, str] = {}
        reason: dict[str, str] = {}
        for f in members:
            m = meta[f]
            calib_level, why = calib[f]
            if m["report_level"]:
                level[f], source[f] = m["report_level"], "families.yaml"
            elif m["scope"] == "annotate_only":
                level[f], source[f] = "none", "families.yaml scope"
            else:
                level[f], source[f] = calib_level, "calibration (report_level not set in families.yaml)"
            if m["scope"] == "annotate_only":
                reason[f] = f"scope: annotate_only; calibration rule: {calib_level} ({why})"
            elif level[f] != calib_level:
                reason[f] = f"families.yaml overrides the calibration rule, which gives {calib_level} ({why})"
                print(f"WARNING: {f}: families.yaml report_level {level[f]}, calibration rule "
                      f"gives {calib_level} ({why})")
            else:
                reason[f] = why

        aggregated = any(level[f] == "group" for f in members)
        if aggregated:
            on = group_choice[0]
            flags = {f: on and f in thresholds for f in members}
        else:
            flags = {f: level[f] == "gene" and gene_choice[f][0] for f in members}

        group_cells = {"f1_all_group_off": f1_cell(group_off, group) if group else "",
                       "f1_all_group_on": (f1_cell(group_on, group) if group_has_on else "NA") if group else ""}
        if aggregated:
            f0 = members[0]
            manifest.append({
                "unit": group, "report_level": "group", "members": ",".join(members),
                "paralog_group": group,
                "cascade_members": ",".join(f for f in members if flags[f]),
                "calibration_level": calib[f0][0], "level_source": source[f0], "reason": reason[f0],
                **chosen((group_on if group_choice[0] else group_off).get(group)),
                "f1_all_gene_off": ",".join(f1_cell(gene_off, f) for f in members),
                "f1_all_gene_on": ",".join(f1_cell(gene_on, f) if f in thresholds else "NA" for f in members),
                **group_cells,
            })
        else:
            for f in members:
                manifest.append({
                    "unit": f, "report_level": level[f], "members": f, "paralog_group": group,
                    "cascade_members": f if flags[f] else "",
                    "calibration_level": calib[f][0], "level_source": source[f], "reason": reason[f],
                    **chosen((gene_on if flags[f] else gene_off).get(f) if level[f] == "gene" else None),
                    "f1_all_gene_off": f1_cell(gene_off, f),
                    "f1_all_gene_on": f1_cell(gene_on, f) if f in thresholds else "NA",
                    **group_cells,
                })

        for f in members:
            if aggregated:
                off_cell = f1_cell(group_off, group)
                on_cell = f1_cell(group_on, group) if group_has_on else "NA"
            else:
                off_cell = f1_cell(gene_off, f)
                on_cell = f1_cell(gene_on, f) if f in thresholds else "NA"
            cascade[f] = {
                "family": f,
                "report_unit": group if aggregated else f,
                "report_level": level[f],
                "diamond_precision": (gene_off.get(f) or {}).get("precision", ""),
                "f1_all_cascade_off": off_cell,
                "f1_all_cascade_on": on_cell,
                "needs_cascade_check": "yes" if flags[f] else "no",
                "hmm_short_read_threshold": thresholds[f] if flags[f] else "",
            }
    return manifest, cascade


MANIFEST_COLUMNS = [
    "unit", "report_level", "members", "paralog_group", "cascade_members",
    "calibration_level", "level_source", "reason",
    "tp", "fp_all", "fn", "precision_all", "recall", "f1_all",
    "f1_all_gene_off", "f1_all_gene_on", "f1_all_group_off", "f1_all_group_on",
]


def write_report_levels(path: Path, rows: list[dict], min_f1: float) -> None:
    header = [
        "osmo_refdb profile report levels: one row per reporting unit (osmotool profile's unit table).",
        "Decided on the calibration reads (positive.calib + negative.train), never on the benchmark (test) reads.",
        f"Rule (min_f1_all = {min_f1}): gene if every member reaches f1_all >= min at gene level; else group "
        "if the paralog group reaches it as one unit; else none.",
        "At each level the cascade (HMM second opinion, osmo_refdb.profile_cascade.tsv) is used only where it "
        "raises f1_all; cascade_members lists where it is on.",
        "f1_all = F1 with every non-background read called as the unit counted as a false positive "
        "(own hard negatives and flanks, other panel genes, other families' hard negatives); "
        "see compute_metrics.py.",
        "report_level is families.yaml's report_level; calibration_level is the rule's. "
        "tp..f1_all are calibration counts at the reported level and cascade (empty for none).",
        "f1_all_gene_*: per member, in members order (NA: no HMM threshold, so no cascade option). "
        "f1_all_group_*: the paralog group as one unit.",
    ]
    with open(path, "w", newline="") as fh:
        for line in header:
            fh.write(f"# {line}\n")
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_report_levels(path: Path) -> list[dict]:
    """Rows of osmo_refdb.profile_report_levels.tsv (comment lines skipped)."""
    with open(path) as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")),
                                   delimiter="\t"))


def summarize_report_units(merged: pd.DataFrame, method: str,
                           units: list[dict]) -> list[dict]:
    """The metrics of summary.tsv per reporting unit of the release's
    osmo_refdb.profile_report_levels.tsv, members scored as one unit. Units
    on none are scored too (so the exclusion can be checked) and keep their
    report_level."""
    n_random_bg = int(merged["_random_bg"].sum())
    rows = []
    for unit in units:
        members = unit["members"].split(",")
        m = _metrics(merged, set(members), n_random_bg)
        if m is None:
            continue
        rows.append({"method": method, "unit": unit["unit"], "report_level": unit["report_level"],
                     "members": unit["members"], "cascade_members": unit["cascade_members"], **m})
    return rows


def verify_unit_counts(cascade_calls: pd.DataFrame, diamond_dir: Path,
                       units: list[dict]) -> None:
    """Check osmotool's <sample>.unit_counts.tsv (osmotool versions that read
    the report-level manifest) against the rebuilt per-read calls: a
    reported unit's raw_count is its members' summed calls (0.5 per paired
    read), an excluded unit must carry no count."""
    files = list(diamond_dir.glob("*.unit_counts.tsv"))
    if not files:
        print("Unit check: no unit_counts.tsv (osmotool without report levels); skipped")
        return
    per_family = (cascade_calls.groupby(["sample", "family"]).size() * 0.5).to_dict()
    members = {u["unit"]: u["members"].split(",") for u in units}
    n_cells = n_bad = 0
    for counts_file in files:
        sample = counts_file.name.replace(".unit_counts.tsv", "")
        df = pd.read_csv(counts_file, sep="\t", comment="#", dtype=str, keep_default_na=False)
        listed = set(df.loc[df["status"] != "not_in_manifest", "unit"])
        if listed != set(members):
            n_bad += 1
            print(f"  MISMATCH {sample}: units {sorted(listed ^ set(members))} "
                  f"not in both the table and the manifest")
        for unit, status, raw in zip(df["unit"], df["status"], df["raw_count"]):
            n_cells += 1
            if status == "excluded":
                ok = raw == "NA"
            elif status == "not_in_manifest":     # e.g. mrpA_mrpB_fused: its own calls
                ok = abs(float(raw) - per_family.get((sample, unit), 0.0)) < 1e-6
            else:
                mine = sum(per_family.get((sample, f), 0.0) for f in members.get(unit, []))
                ok = raw != "NA" and abs(float(raw) - mine) < 1e-6
            if not ok:
                n_bad += 1
                print(f"  MISMATCH {sample} {unit} ({status}): osmotool {raw}")
    print(f"Unit check: osmotool unit_counts.tsv vs rebuilt calls: "
          f"{n_cells - n_bad}/{n_cells} sample x unit cells agree")


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


def summarize_background_fpr(merged: pd.DataFrame, method: str) -> tuple[list[dict], list[dict]]:
    """Background reads -- the shared random bacterial background
    ("random_bg") and an optional real metagenome ("background") -- carry no
    panel gene, so any call on one is a false positive. Returns
    (one row per source, one row per source x called family)."""
    overall, by_family = [], []
    for source in BACKGROUND_SOURCES:
        bg = merged[merged["family_true"] == source]
        if bg.empty:
            continue
        n_reads = len(bg)
        called = bg[bg["family_called"] != "NO_CALL"]
        overall.append({
            "method": method,
            "source": source,
            "n_background_reads": n_reads,
            "n_false_positive_calls": len(called),
            "background_fpr": round(len(called) / n_reads, 6),
            "fp_per_million": round(len(called) / n_reads * 1e6, 1),
        })
        for family, n in called["family_called"].value_counts().items():
            if n == 0:   # categorical dtype lists every category
                continue
            by_family.append({
                "method": method, "source": source, "family_called": family,
                "n_false_positive_calls": int(n), "n_background_reads": n_reads,
                "fp_per_million": round(n / n_reads * 1e6, 1),
            })
    return overall, by_family


def _write_tsv(path: Path, rows: list[dict], label: str) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{label} written to {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=None,
                    help="a directory with reads/, diamond/ and hmm/; or give them with the options below")
    ap.add_argument("--reads", type=Path, default=None, help="default: <results>/reads")
    ap.add_argument("--diamond", type=Path, default=None, help="default: <results>/diamond")
    ap.add_argument("--hmm", type=Path, default=None, help="default: <results>/hmm")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--cascade-out", type=Path, default=None,
                    help="Calibration mode: where to write profile_cascade.tsv "
                         "(osmotool profile's DIAMOND+HMM cascade config). "
                         "Defaults to <out>/profile_cascade.tsv.")
    ap.add_argument("--report-levels-out", type=Path, default=None,
                    help="Calibration mode: where to write profile_report_levels.tsv "
                         "(osmotool profile's unit table). Defaults to "
                         "<out>/profile_report_levels.tsv.")
    ap.add_argument("--min-f1", type=float, default=0.5,
                    help="Calibration mode: the f1_all a unit must reach on the "
                         "calibration reads to be reported at a level (see decide_units).")
    ap.add_argument("--cascade-config", type=Path, default=None,
                    help="Benchmark (test-read) mode: an existing "
                         "profile_cascade.tsv, written earlier from the "
                         "calibration reads. Adds method 'diamond_cascade' "
                         "(osmotool profile's default output) to every "
                         "table; nothing is decided from these reads.")
    ap.add_argument("--report-levels", type=Path, default=None,
                    help="Benchmark mode: the release's profile_report_levels.tsv. "
                         "Adds summary_by_report_unit.tsv and checks osmotool's "
                         "unit_counts.tsv against the rebuilt calls.")
    args = ap.parse_args()

    calibration = args.cascade_config is None
    meta = load_family_meta(args.families)
    families = list(meta)
    groups = {f: m["group"] for f, m in meta.items() if m["group"]}
    report_units = load_report_levels(args.report_levels) if args.report_levels else []
    args.out.mkdir(parents=True, exist_ok=True)
    if args.results is None and not (args.reads and args.diamond and args.hmm):
        ap.error("give --results, or --reads, --diamond and --hmm")
    truth = load_read_truth(args.reads or args.results / "reads")
    diamond_dir = args.diamond or args.results / "diamond"

    diamond_calls = load_diamond_calls(diamond_dir)
    hmm_calls = load_hmm_calls(args.hmm or args.results / "hmm")
    if calibration and not diamond_calls.empty and not any(diamond_dir.glob("*.cascade.hmmsearch.tblout")):
        raise SystemExit(
            f"Calibration mode needs each DIAMOND call's HMM score: run run_benchmark.sh "
            f"with CASCADE=score on these reads (no *.cascade.hmmsearch.tblout in {diamond_dir}).")

    summary_rows = []
    group_rows = []
    best_threshold_rows = []
    length_stratified_rows = []
    background_fpr_rows = []
    background_by_family_rows = []
    all_curve_stats = {}
    merged_by_method: dict[str, pd.DataFrame] = {}

    def add_method(merged: pd.DataFrame, method: str) -> None:
        nonlocal summary_rows, group_rows, length_stratified_rows
        summary_rows += summarize(merged, method, families, groups)
        group_rows += summarize_groups(merged, method, families, groups)
        length_stratified_rows += summarize_by_length(merged, method, families)
        bg_overall, bg_by_family = summarize_background_fpr(merged, method)
        background_fpr_rows.extend(bg_overall)
        background_by_family_rows.extend(bg_by_family)
        merged_by_method[method] = merged

    if not diamond_calls.empty:
        diamond_merged = label_reads(diamond_calls, truth)
        add_method(diamond_merged, "diamond")
        best_threshold_rows += summarize_best_threshold(diamond_merged, "diamond", "bitscore", families)
        all_curve_stats["diamond"] = plot_curves(
            diamond_merged, "diamond", "bitscore", args.out, higher_is_better=True)
        if not calibration:
            cutoffs = load_cascade_config(args.cascade_config)
            cascade_calls = apply_profile_cascade(diamond_calls, diamond_dir, cutoffs)
            verify_cascade_counts(cascade_calls, diamond_dir)
            if report_units:
                verify_unit_counts(cascade_calls, diamond_dir, report_units)
            add_method(label_reads(cascade_calls, truth), "diamond_cascade")
    else:
        print("WARNING: no DIAMOND alignment files found in results/diamond")

    if not hmm_calls.empty:
        hmm_merged = label_reads(hmm_calls, truth)
        add_method(hmm_merged, "hmm")
        best_threshold_rows += summarize_best_threshold(hmm_merged, "hmm", "bitscore", families)
        all_curve_stats["hmm"] = plot_curves(
            hmm_merged, "hmm", "bitscore", args.out, higher_is_better=True)
    else:
        print("WARNING: no HMM tblout files found in results/hmm")

    manifest_rows: list[dict] = []
    cascade_rows: list[dict] = []
    if calibration and not diamond_calls.empty:
        # Every family with an HMM short-read threshold cascaded at once:
        # with the no-cascade tables this scores every option per unit.
        thresholds = {r["family"]: r["best_threshold"]
                      for r in best_threshold_rows if r["method"] == "hmm"}
        all_on = label_reads(apply_profile_cascade(diamond_calls, diamond_dir, thresholds), truth)
        add_method(all_on, "diamond_cascade_all")

        def by(rows, method, key):
            return {r[key]: r for r in rows if r["method"] == method}
        manifest_rows, cascade = decide_units(
            meta, thresholds,
            by(summary_rows, "diamond", "family"), by(summary_rows, "diamond_cascade_all", "family"),
            by(group_rows, "diamond", "paralog_group"), by(group_rows, "diamond_cascade_all", "paralog_group"),
            args.min_f1)
        cascade_rows = list(cascade.values())

    unit_rows = []
    for method, merged in merged_by_method.items():
        if report_units and method != "diamond_cascade_all":
            unit_rows += summarize_report_units(merged, method, report_units)

    _write_tsv(args.out / "summary.tsv", summary_rows, "Per-family precision/recall/F1")
    _write_tsv(args.out / "summary_by_group.tsv", group_rows, "Per-paralog-group precision/recall/F1")
    _write_tsv(args.out / "summary_by_report_unit.tsv", unit_rows,
               "Per-reporting-unit precision/recall/F1 (release report levels)")
    _write_tsv(args.out / "best_threshold_summary.tsv", best_threshold_rows,
               "Per-family best-F1 score-threshold summary")
    _write_tsv(args.out / "summary_by_read_length.tsv", length_stratified_rows,
               "Read-length-stratified precision/recall/F1")
    _write_tsv(args.out / "background_fpr.tsv", background_fpr_rows, "Background false-positive rate")
    _write_tsv(args.out / "background_fpr_by_family.tsv", background_by_family_rows,
               "Background false positives by called family")

    if all_curve_stats:
        print("\nOverall ROC/PR AUC by method:")
        for method, stats in all_curve_stats.items():
            print(f"  {method}: {stats}")

    if not calibration:
        print(f"\nScored the cascade with {args.cascade_config}"
              + (f" and the report levels with {args.report_levels}" if args.report_levels else "")
              + "; nothing decided from these reads.")
        return

    if cascade_rows:
        cascade_path = args.cascade_out or (args.out / "profile_cascade.tsv")
        _write_tsv(cascade_path, cascade_rows, "\nProfile-mode DIAMOND+HMM cascade config")
        flagged = [r["family"] for r in cascade_rows if r["needs_cascade_check"] == "yes"]
        print(f"  {len(flagged)}/{len(cascade_rows)} families cascaded: {', '.join(flagged)}")
    if manifest_rows:
        levels_path = args.report_levels_out or (args.out / "profile_report_levels.tsv")
        write_report_levels(levels_path, manifest_rows, args.min_f1)
        print(f"Profile report levels written to {levels_path}")
        for level in REPORT_LEVELS:
            units = [r["unit"] for r in manifest_rows if r["report_level"] == level]
            print(f"  {level}: {len(units)} units: {', '.join(units)}")


if __name__ == "__main__":
    main()
