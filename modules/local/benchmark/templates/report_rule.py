#!/usr/bin/env python3
"""
report_rule.py: the report-level rule's machinery and the threshold-free
metric, used by both benchmarks (compute_metrics.py for profile mode,
benchmark_annotate.py for annotate mode). The machinery (the gene, group,
none order, the evidence guard, the support test) and the reason vocabulary
are shared, so both manifests read the same way. The settings are not: each
mode has its own threshold, minimum evidence, support and bootstrap
(params.profile_* and params.annotate_*), decided on its own calibration
split, since reads and proteins are different benchmarks for different uses.
Standard library only: the annotate scorer's container has no numpy.

The rule, per paralog group (or per family, for a family without one), on
the calibration split only:
  gene   if every member clears the level test at gene level, with enough
         evidence;
  group  otherwise, if the group as one unit clears it, with enough evidence;
  none   otherwise, member by member, unless that member clears it at gene
         level on its own (then gene).
The level test (passes_level): the point estimate of the metric (f1_all)
reaches min_value AND, with a bootstrap, at least min_support of the
resamples reach it too (level_support). Being unsure pushes a unit down the
order, towards group or none: reporting nothing is better than reporting a
number that does not hold. Without a bootstrap (support None) the test is
the point estimate alone, the rule of releases up to v10.
"Enough evidence" is n_eval_positives >= min_eval_positives: the number of
distinct positive sequences the unit was measured on (summed over a group's
members at group level). A unit below it is none with reason
insufficient_evidence, never below_threshold, whatever its score: a score
measured on a handful of sequences cannot tell a family that works from one
that does not.

The cascade test (cascade_supported, profile mode only): the HMM second
opinion is switched on for a unit only if its point gain in f1_all is
positive and at least min_gain, AND, with a bootstrap, the paired gain is
positive in at least min_support of the resamples (cascade_gain_support).
With min_gain 0 and no bootstrap it is the bare sign of the gain, the rule of
releases up to v10.

The catalog test (unit_test, profile mode only): the level test and the
evidence guard on one reporting unit by itself, on the calibration reads'
catalog calls alone. It decides whether the reads the gene catalog settles
for that unit can be reported, whatever the unit's report level.

The resampling itself (a cluster bootstrap over source sequences) is the
caller's: compute_metrics.py (numpy) for reads, benchmark_annotate.py for
proteins. This module only takes the numbers it produces.

reason, the controlled vocabulary of both manifests:
  meets_threshold        reported at this level by the rule
  below_threshold        measured, and the point estimate fails the rule
  unstable               the point estimate clears the threshold, the
                         bootstrap support does not (reason_detail: both)
  insufficient_evidence  too few positives to decide (min_eval_positives)
  out_of_scope           excluded by families.yaml (scope), not by performance
  override               families.yaml sets another level than the rule's;
                         reason_detail records the rule's level and why
reason_detail is free text: the numbers each choice rests on.
"""

from __future__ import annotations

MEETS = "meets_threshold"
BELOW = "below_threshold"
INSUFFICIENT = "insufficient_evidence"
OUT_OF_SCOPE = "out_of_scope"
OVERRIDE = "override"
UNSTABLE = "unstable"
REASONS = (MEETS, BELOW, UNSTABLE, INSUFFICIENT, OUT_OF_SCOPE, OVERRIDE)


def passes_level(value: float, support: float | None, min_value: float, min_support: float) -> bool:
    """The level test: the point estimate reaches min_value and, when there
    is a bootstrap (support not None), so does at least min_support of the
    resamples."""
    return value >= min_value and (support is None or support >= min_support)


def cascade_supported(gain: float, support: float | None, min_gain: float, min_support: float) -> bool:
    """The cascade test: the point gain (cascade on minus off) is positive and
    at least min_gain and, with a bootstrap, positive in at least
    min_support of the paired resamples. gain is compared as given: pass it
    rounded as the manifest shows it, so the comparison can be redone from
    the file."""
    return gain > 0 and gain >= min_gain and (support is None or support >= min_support)


def _support_text(support: float | None) -> str:
    return "" if support is None else f", support {support:.3f}"


def rule_levels(members: list[str], group: str,
                gene_value: dict[str, float], gene_n: dict[str, int],
                group_value: float | None, group_n: int | None,
                min_value: float, min_n: int, metric: str = "f1_all",
                gene_support: dict[str, float | None] | None = None,
                group_support: float | None = None,
                min_support: float = 0.0) -> dict[str, tuple[str, str, str]]:
    """{family: (level, reason, detail)} for one block: a paralog group's
    members (group set) or a single family (group ""). gene_value/gene_n:
    each member's metric and distinct positive sequences at gene level;
    group_value/group_n: the group as one unit (None without a group). A
    value is a float; a unit with no calls at all counts as 0.

    gene_support / group_support: the fraction of bootstrap resamples in
    which the unit's metric reaches min_value (None, or gene_support None:
    no bootstrap, the point estimate decides alone). A level needs
    passes_level(value, support, min_value, min_support)."""
    gene_support = gene_support or {}

    def enough(n) -> bool:
        return n is not None and n >= min_n

    def ev(n) -> str:
        return f"n_eval_positives {n} < {min_n}"

    def sup(f):
        return gene_support.get(f)

    def need() -> str:
        return f"{metric} >= {min_value}" + (f" in >= {min_support} of resamples" if any(
            sup(f) is not None for f in members) or group_support is not None else "")

    ok = {f: enough(gene_n.get(f)) and passes_level(gene_value[f], sup(f), min_value, min_support)
          for f in members}
    out: dict[str, tuple[str, str, str]] = {}
    if all(ok.values()):
        for f in members:
            out[f] = ("gene", MEETS,
                      f"every member of {group} reaches {metric} {min_value} at gene level "
                      f"({f} {gene_value[f]:.3f}{_support_text(sup(f))})" if group
                      else f"{metric} {gene_value[f]:.3f} >= {min_value}{_support_text(sup(f))}")
        return out

    weak = ", ".join(f"{f} {gene_value[f]:.3f}{_support_text(sup(f))} on {gene_n.get(f) or 0} positives"
                     for f in members if not ok[f])
    if group and enough(group_n) and passes_level(group_value, group_support, min_value, min_support):
        for f in members:
            out[f] = ("group", MEETS, f"{weak} at gene level (needs {need()} on >= {min_n}); "
                                      f"group {metric} {group_value:.3f}{_support_text(group_support)} "
                                      f">= {min_value} on {group_n} positives")
        return out

    tail = ""
    if group:
        if not enough(group_n):
            tail = f"; group {group} {ev(group_n)}"
        elif group_value >= min_value:
            tail = (f"; group {group} {metric} {group_value:.3f} >= {min_value} but support "
                    f"{group_support:.3f} < {min_support}")
        else:
            tail = f"; group {group} {metric} {group_value:.3f} < {min_value}"
    for f in members:
        n, v, s = gene_n.get(f), gene_value[f], sup(f)
        if not enough(n):
            out[f] = ("none", INSUFFICIENT, f"{ev(n)} at gene level ({metric} {v:.3f}){tail}")
        elif ok[f]:
            out[f] = ("gene", MEETS, f"{metric} {v:.3f} >= {min_value}{_support_text(s)} at gene level{tail}")
        elif v >= min_value:
            out[f] = ("none", UNSTABLE, f"{metric} {v:.3f} >= {min_value} but support {s:.3f} < "
                                        f"{min_support} at gene level{tail}")
        else:
            out[f] = ("none", BELOW, f"{metric} {v:.3f} < {min_value}{_support_text(s)} at gene level{tail}")
    return out


def unit_test(value: float, n: int | None, support: float | None, min_value: float, min_n: int,
              min_support: float, metric: str = "f1_all") -> tuple[bool, str, str]:
    """The level test and the evidence guard on one unit by itself, without
    the gene, group, none order: (passes, reason, detail). The profile
    manifest's catalog columns use it, one decision per reporting unit."""
    if n is None or n < min_n:
        return False, INSUFFICIENT, f"n_eval_positives {n or 0} < {min_n} ({metric} {value:.3f})"
    if passes_level(value, support, min_value, min_support):
        return True, MEETS, f"{metric} {value:.3f} >= {min_value}{_support_text(support)}"
    if value >= min_value:
        return False, UNSTABLE, f"{metric} {value:.3f} >= {min_value} but support {support:.3f} < {min_support}"
    return False, BELOW, f"{metric} {value:.3f} < {min_value}{_support_text(support)}"


def rule_header_lines(min_f1: float, min_n: int, n_boot: int, seed: int, min_support: float,
                        unit_word: str = "sequences", min_gain: float | None = None,
                        mode: str = "profile") -> list[str]:
    """The decision rule as manifest header lines, in the same words for both
    manifests, with each mode's own settings and parameter names (mode:
    profile or annotate; min_gain None: annotate, no cascade)."""
    if n_boot:
        level = (f"Level test: f1_all >= {min_f1} on the point estimate AND in >= {min_support} of {n_boot} "
                 f"bootstrap resamples (level_support), resampling source {unit_word} (not reads) within each "
                 f"family's positives and hard negatives, seed {seed}; f1_all_p10/p90: the resamples' 10th and "
                 "90th percentiles. Clearing on the point estimate only gives none for unstable.")
    else:
        level = (f"Level test: f1_all >= {min_f1} on the point estimate alone ({mode}_decision_bootstrap = 0: no "
                 "resampling; level_support, f1_all_p10/p90 and cascade_gain_support are empty).")
    lines = [
        f"Rule: gene if every member passes the level test at gene level; else group if the paralog group "
        f"passes it as one unit; else, member by member, gene if it passes, else none.",
        level,
        f"Minimum evidence ({mode}_min_eval_positives = {min_n}): a level needs that many distinct positive calibration "
        f"{unit_word} (n_eval_positives; a group's summed over its members), else the rule gives none for "
        "insufficient_evidence.",
    ]
    if min_gain is not None:
        lines.append(
            f"Cascade test ({mode}_min_cascade_gain = {min_gain}): the HMM second opinion is on at a unit's level only "
            f"if its point gain in f1_all (cascade_gain, on minus off, 3 decimals) is > 0 and >= {min_gain}"
            + (f", and the paired gain is > 0 in >= {min_support} of the resamples (cascade_gain_support)"
               if n_boot else "")
            + "; the level is judged at the option this picks. cascade_gain NA: no HMM threshold, no option.")
    return lines


def average_precision(scores_pos: list[float], scores_neg: list[float], n_pos: int) -> float:
    """Average precision over a ranking by score, with ties taken together
    (as sklearn's average_precision_score): sum over the distinct scores, from
    the highest, of (recall step) x (precision at that score).

    scores_pos / scores_neg: the scores of the positives and negatives the
    method retrieved, i.e. called as the family. n_pos: every positive,
    retrieved or not. A positive never retrieved contributes nothing, so AP
    is at most the recall at the lowest score. Returns nan without
    positives."""
    if not n_pos:
        return float("nan")
    items = sorted([(s, 1) for s in scores_pos] + [(s, 0) for s in scores_neg], key=lambda x: -x[0])
    ap = tp = fp = 0.0
    prev_recall = 0.0
    i = 0
    while i < len(items):
        s = items[i][0]
        while i < len(items) and items[i][0] == s:
            tp += items[i][1]
            fp += 1 - items[i][1]
            i += 1
        recall = tp / n_pos
        ap += (recall - prev_recall) * tp / (tp + fp)
        prev_recall = recall
    return ap


AP_COLUMN = "ap_vs_own_hard_negatives"
AP_DOC = (f"{AP_COLUMN}: average precision of the family's own score ranking over its positives and its "
          "own hard negatives only (other families' sequences and the background are left out); a sequence "
          "the method did not call as the family is never retrieved and adds nothing, so it is at most "
          "max_recall. ap_random_baseline = n_positive / (n_positive + n_hard_negative), what a random "
          "ranking of the whole set would score.")


def ap_row(scores_pos: list[float], scores_neg: list[float], n_pos: int, n_neg: int) -> dict:
    ap = average_precision(scores_pos, scores_neg, n_pos)
    return {"n_positive": n_pos, "n_hard_negative": n_neg,
            "n_positive_retrieved": len(scores_pos), "n_hard_negative_retrieved": len(scores_neg),
            AP_COLUMN: round(ap, 4) if ap == ap else "",
            "max_recall": round(len(scores_pos) / n_pos, 4) if n_pos else "",
            "ap_random_baseline": round(n_pos / (n_pos + n_neg), 4) if n_pos + n_neg else ""}
