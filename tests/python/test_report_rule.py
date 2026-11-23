"""Unit tests for report_rule.py's decision helpers (standard library only).

Run from the repository root:  python -m unittest discover -s tests/python -v
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "modules/local/benchmark/templates"))
import report_rule as rr  # noqa: E402


class LevelTest(unittest.TestCase):
    def test_point_and_support_clear(self):
        self.assertTrue(rr.passes_level(0.6, 0.95, 0.5, 0.9))

    def test_point_clears_support_does_not(self):
        self.assertFalse(rr.passes_level(0.52, 0.7, 0.5, 0.9))

    def test_point_fails_support_irrelevant(self):
        self.assertFalse(rr.passes_level(0.49, 1.0, 0.5, 0.9))

    def test_boundaries_inclusive(self):
        self.assertTrue(rr.passes_level(0.5, 0.9, 0.5, 0.9))

    def test_no_bootstrap_is_point_estimate(self):
        self.assertTrue(rr.passes_level(0.5, None, 0.5, 0.9))
        self.assertFalse(rr.passes_level(0.499, None, 0.5, 0.9))


class CascadeTest(unittest.TestCase):
    def test_gain_above_minimum_and_supported(self):
        self.assertTrue(rr.cascade_supported(0.039, 0.97, 0.01, 0.9))

    def test_gain_above_minimum_unsupported(self):
        self.assertFalse(rr.cascade_supported(0.027, 0.8, 0.01, 0.9))

    def test_gain_supported_below_minimum(self):
        # opu_permease on v10's calibration reads: +0.006
        self.assertFalse(rr.cascade_supported(0.006, 1.0, 0.01, 0.9))

    def test_gain_at_minimum(self):
        self.assertTrue(rr.cascade_supported(round(0.645 - 0.635, 3), 0.95, 0.01, 0.9))

    def test_bare_sign(self):
        # min_gain 0 and no bootstrap: v10's rule, on for any positive gain, off on a tie
        self.assertTrue(rr.cascade_supported(0.001, None, 0.0, 0.9))
        self.assertFalse(rr.cascade_supported(0.0, None, 0.0, 0.9))
        self.assertFalse(rr.cascade_supported(-0.004, None, 0.0, 0.9))


class RuleLevelsTest(unittest.TestCase):
    members = ["a", "b"]
    n = {"a": 50, "b": 50}

    def levels(self, gv, gs, group_value, group_support):
        return rr.rule_levels(self.members, "g", gv, self.n, group_value, 100, 0.5, 20,
                              gene_support=gs, group_support=group_support, min_support=0.9)

    def test_gene_when_every_member_passes(self):
        out = self.levels({"a": 0.8, "b": 0.7}, {"a": 1.0, "b": 0.95}, 0.75, 1.0)
        self.assertEqual({f: v[:2] for f, v in out.items()}, {"a": ("gene", rr.MEETS), "b": ("gene", rr.MEETS)})

    def test_one_unstable_member_pushes_to_group(self):
        out = self.levels({"a": 0.8, "b": 0.52}, {"a": 1.0, "b": 0.6}, 0.75, 1.0)
        self.assertEqual({f: v[:2] for f, v in out.items()}, {"a": ("group", rr.MEETS), "b": ("group", rr.MEETS)})

    def test_unstable_group_falls_to_members(self):
        out = self.levels({"a": 0.8, "b": 0.52}, {"a": 1.0, "b": 0.6}, 0.55, 0.7)
        self.assertEqual(out["a"][:2], ("gene", rr.MEETS))
        self.assertEqual(out["b"][:2], ("none", rr.UNSTABLE))
        self.assertIn("support 0.600 < 0.9", out["b"][2])
        self.assertIn("0.520 >= 0.5", out["b"][2])

    def test_below_is_not_unstable(self):
        out = self.levels({"a": 0.8, "b": 0.3}, {"a": 1.0, "b": 0.0}, 0.4, 0.0)
        self.assertEqual(out["b"][:2], ("none", rr.BELOW))

    def test_insufficient_evidence_first(self):
        out = rr.rule_levels(["a"], "", {"a": 0.9}, {"a": 5}, None, None, 0.5, 20,
                             gene_support={"a": 1.0}, min_support=0.9)
        self.assertEqual(out["a"][:2], ("none", rr.INSUFFICIENT))

    def test_without_support_matches_point_rule(self):
        out = rr.rule_levels(["a"], "", {"a": 0.52}, {"a": 50}, None, None, 0.5, 20)
        self.assertEqual(out["a"], ("gene", rr.MEETS, "f1_all 0.520 >= 0.5"))

    def test_unstable_in_vocabulary(self):
        self.assertIn(rr.UNSTABLE, rr.REASONS)


if __name__ == "__main__":
    unittest.main()
