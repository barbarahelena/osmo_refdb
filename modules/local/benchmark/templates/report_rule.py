#!/usr/bin/env python3
"""
report_rule.py: the report-level rule and the threshold-free metric shared by
both benchmarks (compute_metrics.py for profile mode, benchmark_annotate.py
for annotate mode), so the two release manifests are decided the same way.
Standard library only: the annotate scorer's container has no numpy.

The rule, per paralog group (or per family, for a family without one), on
the calibration split only:
  gene   if every member clears min_value at gene level, with enough evidence;
  group  otherwise, if the group as one unit clears it, with enough evidence;
  none   otherwise, member by member, unless that member clears at gene level
         on its own (then gene).
"Enough evidence" is n_eval_positives >= min_eval_positives: the number of
distinct positive sequences the unit was measured on (summed over a group's
members at group level). A unit below it is none with reason
insufficient_evidence, never below_threshold, whatever its score: a score
measured on a handful of sequences cannot tell a family that works from one
that does not.

reason, the controlled vocabulary of both manifests:
  meets_threshold        reported at this level by the rule
  below_threshold        measured, and it fails the rule
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
REASONS = (MEETS, BELOW, INSUFFICIENT, OUT_OF_SCOPE, OVERRIDE)


def rule_levels(members: list[str], group: str,
                gene_value: dict[str, float], gene_n: dict[str, int],
                group_value: float | None, group_n: int | None,
                min_value: float, min_n: int, metric: str = "f1_all") -> dict[str, tuple[str, str, str]]:
    """{family: (level, reason, detail)} for one block: a paralog group's
    members (group set) or a single family (group ""). gene_value/gene_n:
    each member's metric and distinct positive sequences at gene level;
    group_value/group_n: the group as one unit (None without a group). A
    value is a float; a unit with no calls at all counts as 0."""
    def enough(n) -> bool:
        return n is not None and n >= min_n

    def ev(n) -> str:
        return f"n_eval_positives {n} < {min_n}"

    ok = {f: enough(gene_n.get(f)) and gene_value[f] >= min_value for f in members}
    out: dict[str, tuple[str, str, str]] = {}
    if all(ok.values()):
        for f in members:
            out[f] = ("gene", MEETS,
                      f"every member of {group} reaches {metric} {min_value} at gene level "
                      f"({f} {gene_value[f]:.3f})" if group else f"{metric} {gene_value[f]:.3f} >= {min_value}")
        return out

    weak = ", ".join(f"{f} {gene_value[f]:.3f} on {gene_n.get(f) or 0} positives"
                     for f in members if not ok[f])
    if group and enough(group_n) and group_value >= min_value:
        for f in members:
            out[f] = ("group", MEETS, f"{weak} at gene level (needs {metric} >= {min_value} on >= {min_n}); "
                                      f"group {metric} {group_value:.3f} >= {min_value} on {group_n} positives")
        return out

    tail = ""
    if group:
        tail = (f"; group {group} {ev(group_n)}" if not enough(group_n)
                else f"; group {group} {metric} {group_value:.3f} < {min_value}")
    for f in members:
        n, v = gene_n.get(f), gene_value[f]
        if not enough(n):
            out[f] = ("none", INSUFFICIENT, f"{ev(n)} at gene level ({metric} {v:.3f}){tail}")
        elif v >= min_value:
            out[f] = ("gene", MEETS, f"{metric} {v:.3f} >= {min_value} at gene level{tail}")
        else:
            out[f] = ("none", BELOW, f"{metric} {v:.3f} < {min_value} at gene level{tail}")
    return out


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
