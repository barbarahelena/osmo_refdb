#!/usr/bin/env python3
"""
check_performance.py -- does a release still perform: its build checks pass
and its accuracy has not dropped against a recorded baseline.

Checks that must hold exactly (any failure fails the run):
  * the release's parts are there (DIAMOND database, HMMs, cascade, report
    levels, markers, catalog, provenance, the benchmark results);
  * every family has an HMM with a gathering cutoff, a reporting unit in both
    report-level manifests (profile and annotate), each with a reason from
    the fixed vocabulary, and a row in the benchmark; every marker is in
    markers/;
  * the per-read rebuilds agree with osmotool's own counts, in the read
    benchmark (results/metrics/verification.tsv) and the catalog benchmark
    (results/catalog_benchmark/verification.tsv), and annotate's hmm best hits
    with its calls (results/annotate/verification.tsv);
  * the catalog does not make things worse: median f1_all over the reported
    units with every sequence catalogued (A) >= with none of the scored ones (B).

Accuracy against the baseline (--baseline, a TSV of metric and value): each
metric of the read benchmark (diamond_cascade and hmm, per reporting unit:
f1_all, precision_all, recall), the annotate benchmark (per family: f1,
precision, recall), the catalog benchmark (catalog_A and catalog_B, per unit)
and the real-genome reads (diamond, hmm, catalog, per target family: recall,
precision). A metric more than --tolerance below its baseline is a regression;
more than --tolerance above, an improvement to record. Regressions fail the run
unless the release was built from another UniProt release than the baseline
(its references differ, so its numbers can): then they are warnings, unless
--strict.

Without a baseline it writes the run's metrics as one, to --new-baseline if
given (the nf-test test: its output folder, so a test leaves nothing in the
clone), else to --baseline; --write-baseline overwrites --baseline. Review it
before committing it. Every run writes its metrics to --out.

Standard library only.

Usage:
  python tests/check_performance.py --release-dir <release folder> \\
      --families assets/test/families.yaml --markers assets/test/markers.yaml \\
      --baseline tests/baseline.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from pathlib import Path

REQUIRED = [
    "osmo_refdb.dmnd", "osmo_refdb.diamond_cutoffs.tsv", "osmo_refdb.train_refs.faa",
    "osmo_refdb.profile_cascade.tsv", "osmo_refdb.profile_report_levels.tsv",
    "osmo_refdb.annotate_report_levels.tsv",
    "osmo_refdb.profile_excluded_families.txt", "osmo_refdb.annotate_excluded_families.txt",
    "hmms/osmo_refdb.hmm", "hmms/cutoff_manifest.tsv", "refs/manifest.tsv",
    "markers/osmo_markers.dmnd", "markers/osmo_markers.manifest.tsv", "markers/osmo_markers.train_refs.faa",
    "PROVENANCE.json", "software_versions.yml",
    "results/metrics/summary_by_report_unit.tsv", "results/metrics/verification.tsv",
    "results/metrics/average_precision.tsv",
    "results/annotate/summary.tsv", "results/annotate/summary_by_group.tsv",
    "results/annotate/average_precision.tsv",
    "results/annotate_calib/summary.tsv", "results/annotate_calib/summary_by_group.tsv",
    "results/calib_metrics/decision_sensitivity.tsv",
    "results/catalog_benchmark/summary_by_report_unit.tsv", "results/catalog_benchmark/verification.tsv",
]
HIGHER_IS_BETTER = ("f1", "f1_all", "precision", "precision_all", "recall")
REPORT_LEVEL_FILES = ("osmo_refdb.profile_report_levels.tsv", "osmo_refdb.annotate_report_levels.tsv")
REASONS = ("meets_threshold", "below_threshold", "unstable", "insufficient_evidence", "out_of_scope",
           "override")
# the decision bootstrap, in both manifests (cascade_gain*: NA in annotate's, which has no cascade)
ROBUST_COLUMNS = ("f1_all_p10", "f1_all_p90", "level_support", "cascade_gain", "cascade_gain_support")
# the catalog decision, profile manifest only: a build always runs it, so NA is a failure here
CATALOG_LEVEL_FILE = "osmo_refdb.profile_report_levels.tsv"
CATALOG_REASONS = ("meets_threshold", "below_threshold", "unstable", "insufficient_evidence", "out_of_scope")
DECISION_PARAMS = ("profile_min_report_f1", "profile_min_eval_positives", "profile_min_cascade_gain",
                   "profile_decision_support", "profile_decision_bootstrap",
                   "annotate_min_report_f1", "annotate_min_eval_positives", "annotate_decision_support",
                   "annotate_decision_bootstrap", "decision_bootstrap_seed")


def yaml_names(path: Path, section: str) -> list[str]:
    """The `- name:` entries of a top-level YAML section (no YAML library)."""
    names, inside = [], False
    for line in path.read_text().splitlines():
        if re.match(r"^[A-Za-z_]+:", line):
            inside = line.startswith(f"{section}:")
            continue
        m = re.match(r"^\s*-\s+name:\s*['\"]?([^'\"\s#]+)", line)
        if inside and m:
            names.append(m.group(1))
    return names


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))


def num(v) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def hmm_cutoffs(path: Path) -> dict[str, bool]:
    """{model name: has a GA line} of a (concatenated) HMMER3 file."""
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith("NAME "):
            name = line.split()[1]
            out[name] = False
        elif line.startswith("GA ") and name:
            out[name] = True
    return out


class Checks:
    def __init__(self):
        self.failures: list[str] = []
        self.passed = 0

    def require(self, ok: bool, what: str) -> bool:
        if ok:
            self.passed += 1
        else:
            self.failures.append(what)
        return ok


def build_checks(rel: Path, families: list[str], markers: list[str], c: Checks) -> None:
    for f in REQUIRED:
        c.require((rel / f).exists(), f"missing {f}")
    c.require(any((rel / "catalogs").glob("*/catalog.json")), "no catalog under catalogs/")
    # a shipped catalog names the release it was built for, as PROVENANCE.json does (not the
    # directory Nextflow staged the release under)
    prov = rel / "PROVENANCE.json"
    if prov.exists():
        name = json.loads(prov.read_text()).get("release")
        for cj in sorted((rel / "catalogs").glob("*/catalog.json")):
            got = json.loads(cj.read_text()).get("release")
            c.require(got == name, f"{cj.relative_to(rel)}: release is {got!r}, PROVENANCE.json says {name!r}")

    if (rel / "hmms/osmo_refdb.hmm").exists():
        ga = hmm_cutoffs(rel / "hmms/osmo_refdb.hmm")
        for fam in families:
            c.require(ga.get(fam, False), f"{fam}: no HMM with a GA cutoff in hmms/osmo_refdb.hmm")
    for levels in REPORT_LEVEL_FILES:
        if not (rel / levels).exists():
            continue
        members = set()
        rows = read_tsv(rel / levels)
        for col in ROBUST_COLUMNS:
            c.require(bool(rows) and col in rows[0], f"{levels}: no column {col}")
        for r in rows:
            members.update(m for m in r.get("members", "").split(",") if m)
            c.require(r.get("reason") in REASONS,
                      f"{levels}: unit {r.get('unit')}: reason {r.get('reason')!r} is not one of {', '.join(REASONS)}")
            if levels == CATALOG_LEVEL_FILE:
                scoped = r.get("reason") == "out_of_scope"
                c.require(r.get("catalog_reported") == "NA" if scoped else r.get("catalog_reported") in ("yes", "no"),
                          f"{levels}: unit {r.get('unit')}: catalog_reported {r.get('catalog_reported')!r} "
                          + ("is not NA for a unit out of scope" if scoped else
                             "is not yes or no (no catalog run on the calibration reads?)"))
                c.require(r.get("catalog_reason") in CATALOG_REASONS,
                          f"{levels}: unit {r.get('unit')}: catalog_reason {r.get('catalog_reason')!r}")
        for fam in families:
            c.require(fam in members, f"{fam}: in no reporting unit ({levels})")
    if (rel / "markers/osmo_markers.manifest.tsv").exists():
        listed = {r["marker"] for r in read_tsv(rel / "markers/osmo_markers.manifest.tsv")}
        for m in markers:
            c.require(m in listed, f"marker {m}: not in markers/osmo_markers.manifest.tsv")
    if (rel / "PROVENANCE.json").exists():
        try:
            prov = json.loads((rel / "PROVENANCE.json").read_text())
            for key in ("inputs", "outputs"):
                c.require(key in prov, f"PROVENANCE.json has no {key}")
            rule = prov.get("inputs", {}).get("decision_rule", {})
            for key in DECISION_PARAMS:
                c.require(rule.get(key, "not recorded") != "not recorded",
                          f"PROVENANCE.json inputs.decision_rule does not record {key}")
        except json.JSONDecodeError as e:
            c.require(False, f"PROVENANCE.json does not parse: {e}")

    for v in ("results/metrics/verification.tsv", "results/catalog_benchmark/verification.tsv"):
        if (rel / v).exists():
            rows = read_tsv(rel / v)
            c.require(bool(rows), f"{v} is empty")
            for r in rows:
                what = r.get("check") or r.get("sample") or "?"
                c.require(r.get("cells_disagree", r.get("disagree")) == "0",
                          f"{v}: {what} disagrees with osmotool ({r})")

    v = "results/annotate/verification.tsv"
    if (rel / v).exists():
        for r in read_tsv(rel / v):
            if r.get("method") == "hmm":      # diamond's best hit is taken before its cutoffs: not exact
                c.require(r.get("proteins_disagree") == "0",
                          f"{v}: hmm best hits disagree with osmotool's calls ({r})")

    if (rel / "results/metrics/summary.tsv").exists():
        seen = {r["family"] for r in read_tsv(rel / "results/metrics/summary.tsv")}
        for fam in families:
            c.require(fam in seen, f"{fam}: no row in results/metrics/summary.tsv")

    cat = rel / "results/catalog_benchmark/summary_by_report_unit.tsv"
    if cat.exists():
        rows = read_tsv(cat)

        def med(method):
            vals = [num(r["f1_all"]) or 0.0 for r in rows
                    if r["method"] == method and r.get("report_level") != "none"]
            return statistics.median(vals) if vals else None
        a, b = med("catalog_A"), med("catalog_B")
        if a is not None and b is not None:
            c.require(a >= b, f"catalog benchmark: median f1_all with every sequence catalogued ({a:.3f}) "
                              f"below that with none catalogued ({b:.3f})")


def metrics(rel: Path, families: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}

    def add(key, value):
        v = num(value)
        if v is not None:
            out[key] = v
    p = rel / "results/metrics/summary_by_report_unit.tsv"
    if p.exists():
        for r in read_tsv(p):
            if r["method"] in ("diamond_cascade", "hmm"):
                for m in ("f1_all", "precision_all", "recall"):
                    add(f"read.{r['method']}.{r['unit']}.{m}", r.get(m))
    p = rel / "results/annotate/summary.tsv"
    if p.exists():
        for r in read_tsv(p):
            for m in ("f1", "precision", "recall"):
                add(f"annotate.{r['method']}.{r['family']}.{m}", r.get(m))
    p = rel / "results/catalog_benchmark/summary_by_report_unit.tsv"
    if p.exists():
        for r in read_tsv(p):
            if r["method"] in ("catalog_A", "catalog_B"):
                for m in ("f1_all", "precision_all", "recall"):
                    add(f"catalog.{r['method']}.{r['unit']}.{m}", r.get(m))
    for p in sorted((rel / "real_genome_validation").glob("*.real_recall_comparison.tsv")):
        genome = p.name.split(".")[0]
        for r in read_tsv(p):
            if r["family"] in families:
                for m in ("recall", "precision"):
                    add(f"real.{genome}.{r['method']}.{r['family']}.{m}", r.get(m))
    return out


def uniprot_release(rel: Path) -> str:
    p = rel / "refs/manifest.tsv"
    if not p.exists():
        return "unknown"
    vals = {r.get("uniprot_release", "") for r in read_tsv(p)} - {""}
    return ",".join(sorted(vals)) or "unknown"


def write_metrics(path: Path, values: dict[str, float], release: str) -> None:
    with open(path, "w") as fh:
        fh.write(f"# uniprot_release\t{release}\n")
        fh.write("metric\tvalue\n")
        for k in sorted(values):
            fh.write(f"{k}\t{values[k]:.4f}\n")


def read_baseline(path: Path) -> tuple[dict[str, float], str]:
    release = "unknown"
    with open(path) as fh:
        for line in fh:
            if line.startswith("# uniprot_release\t"):
                release = line.rstrip("\n").split("\t", 1)[1]
    return {r["metric"]: float(r["value"]) for r in read_tsv(path)}, release


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release-dir", type=Path, required=True)
    ap.add_argument("--families", type=Path, required=True)
    ap.add_argument("--markers", type=Path, required=True)
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--tolerance", type=float, default=0.05, help="absolute, on every metric (0-1)")
    ap.add_argument("--out", type=Path, default=None, help="this run's metrics (default: <release>/performance.tsv)")
    ap.add_argument("--write-baseline", action="store_true", help="(over)write --baseline with this run's metrics")
    ap.add_argument("--new-baseline", type=Path, default=None,
                    help="without a baseline, write this run's metrics here instead of to --baseline "
                         "(review it, then copy it to --baseline)")
    ap.add_argument("--strict", action="store_true", help="regressions fail even across UniProt releases")
    args = ap.parse_args()

    rel = args.release_dir
    families = yaml_names(args.families, "families")
    markers = yaml_names(args.markers, "families")
    print(f"release {rel}: {len(families)} families, {len(markers)} markers")

    c = Checks()
    build_checks(rel, families, markers, c)
    values = metrics(rel, families)
    release = uniprot_release(rel)
    out = args.out or rel / "performance.tsv"
    write_metrics(out, values, release)
    print(f"build checks: {c.passed} passed, {len(c.failures)} failed; {len(values)} metrics -> {out}")
    for f in c.failures:
        print(f"  FAIL {f}")
    sysrep = rel / "real_genome_validation/reference_systems_check.txt"
    if sysrep.exists():
        last = sysrep.read_text().strip().splitlines()[-1:]
        print(f"reference systems (informative, not checked): {last[0] if last else ''}")

    regressions, warnings = [], []
    if args.write_baseline or not args.baseline.exists():
        target = args.baseline if args.write_baseline or args.new_baseline is None else args.new_baseline
        write_metrics(target, values, release)
        print(f"no baseline compared: wrote {target} from this run "
              f"(review it, then commit it as {args.baseline})")
    else:
        base, base_release = read_baseline(args.baseline)
        other_release = base_release != release
        improved, missing = [], []
        for k, b in sorted(base.items()):
            v = values.get(k)
            if v is None:
                missing.append(k)
            elif v < b - args.tolerance:
                (warnings if other_release and not args.strict else regressions).append(f"{k}: {v:.3f} (baseline {b:.3f})")
            elif v > b + args.tolerance:
                improved.append(f"{k}: {v:.3f} (baseline {b:.3f})")
        new = sorted(set(values) - set(base))
        print(f"against {args.baseline} (UniProt {base_release}; this run {release}), tolerance {args.tolerance}: "
              f"{len(base) - len(missing)} compared, {len(regressions)} regressions, {len(warnings)} warnings, "
              f"{len(improved)} improvements, {len(missing)} missing, {len(new)} new")
        for r in regressions:
            print(f"  REGRESSION {r}")
        for w in warnings:
            print(f"  WARNING (other UniProt release) {w}")
        for i in improved:
            print(f"  improved {i}")
        for m in missing:
            regressions.append(f"{m}: missing from this run")
            print(f"  MISSING {m}")
        if improved or new:
            print("  (record them with --write-baseline once reviewed)")

    if c.failures or regressions:
        print("RESULT: FAIL")
        return 1
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
