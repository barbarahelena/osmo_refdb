"""Tests for the bootstrapped report-level and cascade decision.

compute_metrics.py needs numpy, pandas, matplotlib, scikit-learn and pyyaml
(the metrics container); benchmark_annotate.py needs pyyaml. Run from the
repository root:  python -m unittest discover -s tests/python -v

  - profile_decision_bootstrap 0 with profile_min_cascade_gain 0 reproduces v10's manifests
    from v10's calibration tables, cell for cell (data/v10_calibration);
  - the bootstrap's clusters add up to the summary tables it decides on;
  - the same seed gives the same manifest, byte for byte;
  - the cascade comparison is paired: an option identical to off is never
    supported;
  - the annotate protein bootstrap: the same checks.
"""
import csv
import filecmp
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = Path(__file__).resolve().parent / "data" / "v10_calibration"
sys.path.insert(0, str(REPO / "modules/local/benchmark/metrics/templates"))
sys.path.insert(0, str(REPO / "modules/local/benchmark/templates"))
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp())

try:
    import numpy as np
    import pandas as pd
    import compute_metrics as cm
except ImportError as e:          # pragma: no cover
    cm = None
    MISSING = str(e)
import benchmark_annotate as ba  # noqa: E402


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))


def ints(row):
    for k in ("tp", "fn", "fp_all"):
        if row.get(k, "") != "":
            row[k] = int(row[k])
    return row


@unittest.skipIf(cm is None, "compute_metrics.py's packages are not installed")
class V10Replay(unittest.TestCase):
    """The audit: without the bootstrap and the minimum gain, the decision is v10's."""

    def decide(self, min_gain):
        summary = [ints(r) for r in read_tsv(DATA / "summary.tsv")]
        groups = [ints(r) for r in read_tsv(DATA / "summary_by_group.tsv")]
        thresholds = {r["family"]: float(r["best_threshold"])
                      for r in read_tsv(DATA / "best_threshold_summary.tsv") if r["method"] == "hmm"}

        def by(rows, method, key):
            return {r[key]: r for r in rows if r["method"] == method}
        meta = cm.load_family_meta(DATA / "families.yaml")
        return cm.decide_units(meta, thresholds,
                               by(summary, "diamond", "family"), by(summary, "diamond_cascade_all", "family"),
                               by(groups, "diamond", "paralog_group"),
                               by(groups, "diamond_cascade_all", "paralog_group"),
                               0.5, 0, None, min_gain, 0.9)

    def test_reproduces_v10(self):
        manifest, cascade, _ = self.decide(0.0)
        old = read_tsv(DATA / "osmo_refdb.profile_report_levels.tsv")
        self.assertEqual([r["unit"] for r in manifest], [r["unit"] for r in old])
        for o, n in zip(old, manifest):
            for k in o:
                if k == "reason":          # v10's free text, now reason_detail; reason is the code
                    continue
                self.assertEqual(str(n[k]), o[k], f"{o['unit']} {k}")
        for o in read_tsv(DATA / "osmo_refdb.profile_cascade.tsv"):
            for k in o:
                self.assertEqual(str(cascade[o["family"]][k]), o[k], f"{o['family']} {k}")

    def test_min_cascade_gain_on_v10(self):
        _, cascade, sensitivity = self.decide(0.01)
        switched = sorted(r["unit"] for r in sensitivity if r["cascade_differs"] == "yes")
        self.assertEqual(switched, ["galE", "gshA", "opu_permease"])
        self.assertEqual(cascade["opuAB"]["needs_cascade_check"], "no")
        self.assertEqual(cascade["proX"]["needs_cascade_check"], "yes")


def synthetic_reads(seed=7):
    """Calibration reads of three families (a, b: paralog group g; c alone)
    plus background, and their calls without and with the cascade."""
    rng = np.random.default_rng(seed)
    truth, off, on = [], [], []
    rid = 0
    for fam, sib in (("a", "b"), ("b", "a"), ("c", None)):
        for kind, n_seq, reads in (("positive", 40, 15), ("negative", 40, 8)):
            for s in range(n_seq):
                uid = f"{fam}{kind[0]}{s}"
                quality = rng.uniform(0.4, 0.95)             # sequences differ: reads of one are correlated
                for _ in range(reads):
                    rid += 1
                    genuine = kind == "positive" and rng.random() > 0.1
                    truth.append({"read_id": f"r{rid}/1", "family": fam, "label": kind, "uniprot_id": uid,
                                  "overlaps_cds": int(genuine), "read_length": 150})
                    x = rng.random()
                    if genuine:
                        call = fam if x < quality else (sib if sib and x < quality + 0.15 else None)
                    else:
                        call = fam if x < 0.12 else None
                    if call:
                        off.append({"sample": "s1", "read_id": f"r{rid}/1", "family": call,
                                    "bitscore": float(rng.uniform(30, 90))})
                        # the cascade drops some calls, more of them false
                        if rng.random() > (0.05 if genuine else 0.5):
                            on.append(off[-1])
    for i in range(200):
        rid += 1
        truth.append({"read_id": f"r{rid}/1", "family": "random_bg", "label": "negative", "uniprot_id": f"bg{i}",
                      "overlaps_cds": 0, "read_length": 150})
        if i % 20 == 0:
            off.append({"sample": "s1", "read_id": f"r{rid}/1", "family": "c", "bitscore": 31.0})
    truth = pd.DataFrame(truth)
    return (cm.label_reads(pd.DataFrame(off), truth), cm.label_reads(pd.DataFrame(on), truth))


def synthetic_catalog(merged, seed=9):
    """Catalog-only calls on the same reads, as BENCHMARK_CATALOGMAP's
    catalog_hits.tsv would give them: c's settled reads are called right, a
    and b are called on their own hard negatives as often as on their
    positives. Returns (catalog_hits rows, the truth table)."""
    rng = np.random.default_rng(seed)
    truth = merged[["read_id", "family_true", "label", "uniprot_id", "overlaps_cds", "read_length"]].rename(
        columns={"family_true": "family"}).astype({"family": str})
    hits = []
    for r in truth.itertuples():
        if r.family == "random_bg":
            continue
        genuine = r.label == "positive" and r.overlaps_cds
        if r.family == "c":
            call = genuine and rng.random() < 0.9
        else:
            call = rng.random() < (0.5 if genuine or r.label == "negative" else 0.0)
        read, mate = r.read_id.rsplit("/", 1)
        if call:
            hits.append({"read": read, "mate": mate, "decision": "call", "label": r.family})
        elif rng.random() < 0.1:
            hits.append({"read": read, "mate": mate, "decision": "ambiguous", "label": "a,b"})
    return hits, truth


@unittest.skipIf(cm is None, "compute_metrics.py's packages are not installed")
class ProfileBootstrap(unittest.TestCase):
    families = ["a", "b", "c"]
    groups = {"a": "g", "b": "g"}
    meta = {"a": {"group": "g"}, "b": {"group": "g"}, "c": {"group": ""}}
    for _m in meta.values():
        _m.update(report_level="", report_level_reason="", scope="")
    thresholds = {"a": 40.0, "b": 40.0, "c": 40.0}
    units = {("gene", "a"): {"a"}, ("gene", "b"): {"b"}, ("gene", "c"): {"c"}, ("group", "g"): {"a", "b"}}

    @classmethod
    def setUpClass(cls):
        cls.off, cls.on = synthetic_reads()

    def tables(self):
        def rows(merged, fn):
            return fn(merged, "x", self.families, self.groups)
        return ({r["family"]: r for r in rows(self.off, cm.summarize)},
                {r["family"]: r for r in rows(self.on, cm.summarize)},
                {r["paralog_group"]: r for r in rows(self.off, cm.summarize_groups)},
                {r["paralog_group"]: r for r in rows(self.on, cm.summarize_groups)})

    def manifest(self, seed, n_boot=200, path=None):
        boot, point = cm.cluster_bootstrap({"off": self.off, "on": self.on}, self.units, n_boot, seed)
        go, gn, ro, rn = self.tables()
        m, c, s = cm.decide_units(self.meta, self.thresholds, go, gn, ro, rn, 0.5, 20, boot, 0.01, 0.9)
        if path:
            cm.write_report_levels(path, m, 0.5, 20, n_boot, seed, 0.9, 0.01)
        return m, c, s, point

    def test_clusters_add_up_to_summary(self):
        go, gn, ro, rn = self.tables()
        _, point = cm.cluster_bootstrap({"off": self.off, "on": self.on}, self.units, 10, 1)
        cm.check_bootstrap_points(point, {
            (("gene", f), "off"): go.get(f) for f in self.families} | {
            (("gene", f), "on"): gn.get(f) for f in self.families} | {
            (("group", "g"), "off"): ro["g"], (("group", "g"), "on"): rn["g"]})

    def test_same_seed_same_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            self.manifest(42, path=Path(d) / "one.tsv")
            self.manifest(42, path=Path(d) / "two.tsv")
            self.manifest(43, path=Path(d) / "other.tsv")
            self.assertTrue(filecmp.cmp(Path(d) / "one.tsv", Path(d) / "two.tsv", shallow=False))
            self.assertFalse(filecmp.cmp(Path(d) / "one.tsv", Path(d) / "other.tsv", shallow=False))

    def test_columns_filled(self):
        manifest, cascade, sensitivity, _ = self.manifest(42)
        for row in manifest:
            self.assertIn(row["reason"], cm.report_rule.REASONS)
            self.assertNotEqual(row["level_support"], "")
            self.assertLessEqual(row["f1_all_p10"], row["f1_all_p90"])
        self.assertEqual({r["unit"] for r in sensitivity}, {r["unit"] for r in manifest})
        self.assertTrue(all(c["cascade_gain"] != "" for c in cascade.values()))

    def catalog(self):
        """The catalog calls, read back from a catalog_hits.tsv as the pipeline does."""
        hits, truth = synthetic_catalog(self.off)
        with tempfile.TemporaryDirectory() as d:
            pd.DataFrame(hits).to_csv(Path(d) / "s1.catalog_hits.tsv", sep="\t", index=False)
            settled = cm.catalog_decisions(Path(d), ["s1"])
        self.assertEqual(len(settled), len(hits))
        return cm.label_reads(cm.catalog_only_calls(settled), truth)

    def test_catalog_leaves_level_and_cascade(self):
        """The catalog is one more configuration of the same draws: the off/on resamples, and so
        every level and cascade, are what they are without it."""
        cat = self.catalog()
        plain, _ = cm.cluster_bootstrap({"off": self.off, "on": self.on}, self.units, 100, 42)
        both, _ = cm.cluster_bootstrap({"off": self.off, "on": self.on, "catalog": cat}, self.units, 100, 42)
        for key in self.units:
            for c in ("off", "on"):
                self.assertTrue(np.array_equal(plain[key][c], both[key][c]), (key, c))
        go, gn, ro, rn = self.tables()
        gc = {r["family"]: r for r in cm.summarize(cat, "x", self.families, self.groups)}
        rc = {r["paralog_group"]: r for r in cm.summarize_groups(cat, "x", self.families, self.groups)}
        without, cascade0, _ = cm.decide_units(self.meta, self.thresholds, go, gn, ro, rn, 0.5, 20, plain, 0.01, 0.9)
        with_cat, cascade1, sens = cm.decide_units(self.meta, self.thresholds, go, gn, ro, rn, 0.5, 20, both,
                                                   0.01, 0.9, gc, rc)
        self.assertEqual(cascade0, cascade1)
        for a, b in zip(without, with_cat):
            self.assertEqual({k: v for k, v in b.items() if k not in cm.CATALOG_COLUMNS},
                             {k: v for k, v in a.items() if k not in cm.CATALOG_COLUMNS})
            self.assertEqual(a["catalog_reported"], "NA")
        by_unit = {r["unit"]: r for r in with_cat}
        self.assertEqual(by_unit["c"]["catalog_reported"], "yes")
        for row in with_cat:
            self.assertIn(row["catalog_reported"], ("yes", "no"))
            self.assertIn(row["catalog_reason"], cm.report_rule.REASONS)
            self.assertNotEqual(row["catalog_level_support"], "")
        self.assertTrue(all(r["catalog_differs"] in ("yes", "no") for r in sens))
        with tempfile.TemporaryDirectory() as d:
            cm.write_report_levels(Path(d) / "m.tsv", with_cat, 0.5, 20, 100, 42, 0.9, 0.01)
            back = cm.load_report_levels(Path(d) / "m.tsv")
        self.assertEqual([r["catalog_reported"] for r in back], [r["catalog_reported"] for r in with_cat])

    def test_paired_comparison(self):
        """An option identical to off has gain 0 in every resample: never supported."""
        boot, _ = cm.cluster_bootstrap({"off": self.off, "on": self.off}, self.units, 50, 3)
        for key in self.units:
            self.assertTrue(np.array_equal(boot[key]["off"], boot[key]["on"]))

    def test_resamples_are_by_sequence(self):
        """Resampling sequences, not reads, spreads f1_all wider than a read-level resample would:
        here the reads of one sequence share its quality."""
        boot, _ = cm.cluster_bootstrap({"off": self.off}, {("gene", "c"): {"c"}}, 400, 5)
        cluster_sd = boot[("gene", "c")]["off"].std()
        own = self.off[self.off["family_true"] == "c"]
        reads = own[["is_genuine_positive", "family_called"]].astype(str).to_numpy()
        rng = np.random.default_rng(5)
        f1 = []
        for _ in range(400):
            r = reads[rng.integers(0, len(reads), len(reads))]
            genuine, called = r[:, 0] == "True", r[:, 1] == "c"
            tp, fn, fp = (genuine & called).sum(), (genuine & ~called).sum(), (~genuine & called).sum()
            f1.append(2 * tp / (2 * tp + fn + fp))
        self.assertGreater(cluster_sd, np.std(f1))


class AnnotateBootstrap(unittest.TestCase):
    families = ["a", "b", "c"]
    groups = {"g": {"a", "b"}}

    def proteins(self):
        import random
        rnd = random.Random(11)
        rows, calls = [], {}
        for fam, sib in (("a", "b"), ("b", "a"), ("c", None)):
            for kind in ("positive", "negative"):
                for i in range(60):
                    pid = f"{fam}{kind[0]}{i}"
                    rows.append({"kind": kind, "source_family": fam, "protein_id": pid})
                    x = rnd.random()
                    if kind == "positive":
                        calls[pid] = fam if x < 0.6 else (sib if sib and x < 0.75 else None)
                    elif x < 0.1:
                        calls[pid] = fam
        for i in range(50):
            rows.append({"kind": "background", "source_family": "random_bg", "protein_id": f"bg{i}"})
            if i % 10 == 0:
                calls[f"bg{i}"] = "c"
        return rows, {k: v for k, v in calls.items() if v}

    def units(self):
        return {**{("gene", f): {f} for f in self.families}, ("group", "g"): self.groups["g"]}

    def test_point_counts_match_score(self):
        rows, calls = self.proteins()
        _, point = ba.protein_bootstrap(rows, calls, self.units(), 5, 1)
        for key, members in self.units().items():
            s = ba.score(rows, calls, members, 50)
            self.assertEqual(point[key], (s["tp"], s["fn"], s["fp_all"]), key)

    def test_same_seed_same_resamples(self):
        rows, calls = self.proteins()
        one, _ = ba.protein_bootstrap(rows, calls, self.units(), 100, 42)
        two, _ = ba.protein_bootstrap(rows, calls, self.units(), 100, 42)
        other, _ = ba.protein_bootstrap(rows, calls, self.units(), 100, 43)
        self.assertEqual(one, two)
        self.assertNotEqual(one, other)

    def test_levels_with_and_without_bootstrap(self):
        rows, calls = self.proteins()
        fam_cfg = [{"name": f, **({"paralog_group": "g"} if f in "ab" else {})} for f in self.families]
        gene = {f: {"family": f, **ba.score(rows, calls, {f}, 50)} for f in self.families}
        group = {"g": {"group": "g", **ba.score(rows, calls, self.groups["g"], 50)}}
        boot, _ = ba.protein_bootstrap(rows, calls, self.units(), 300, 42)
        point = ba.decide_levels(fam_cfg, gene, group, "hmm", 0.5, 20)
        robust = ba.decide_levels(fam_cfg, gene, group, "hmm", 0.5, 20, boot=boot, min_support=0.9)
        for row in point:
            self.assertEqual(row["level_support"], "")
        for row in robust:
            self.assertIn(row["reason"], ba.report_rule.REASONS)
            self.assertNotEqual(row["level_support"], "")

    def test_percentile_matches_numpy(self):
        if cm is None:
            self.skipTest("numpy not installed")
        xs = list(np.random.default_rng(0).random(101))
        for q in (0, 10, 50, 90, 100):
            self.assertAlmostEqual(ba.percentile(xs, q), float(np.percentile(xs, q)))


if __name__ == "__main__":
    unittest.main()
