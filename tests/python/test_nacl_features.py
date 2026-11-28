"""Tests for the NaCl model's features restricted to annotate's reporting units
(nacl_model/features.py, stdlib only). Run from the repository root:
  python -m unittest discover -s tests/python -v

  - read_report_levels on v11's annotate_report_levels.tsv: 30 families at
    level gene, opu_permease as one group, 11 families not reported;
  - a group feature is the sum of its members' raw_count, the members are not
    features themselves, and the breadth summaries count them;
  - without groups the features are as before (schema 1);
  - feature_names puts the groups after the families.
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LEVELS = Path(__file__).resolve().parent / "data" / "v11_annotate" / "osmo_refdb.annotate_report_levels.tsv"
sys.path.insert(0, str(REPO / "modules/local/nacl/model/templates/nacl_model"))
import features as F  # noqa: E402

NOT_REPORTED_V11 = {"opuAA", "opuAC", "opuBA", "opuBC", "opuCA", "opuCC",
                    "mrpC", "mrpD", "trkH", "ktrB", "ktrD"}


def write_gene_counts(prefix, counts):
    lines = ["# total_proteins\t4000", "# gc\t0.5", "# n_contigs\t3", "gene\traw_count"]
    lines += [f"{g}\t{c}" for g, c in counts.items()] + ["betL_decoy\t5"]
    Path(f"{prefix}.gene_counts.tsv").write_text("\n".join(lines) + "\n")


class ReportLevels(unittest.TestCase):
    def test_v11_units(self):
        families, groups, excluded = F.read_report_levels(LEVELS)
        self.assertEqual(len(families), 30)
        self.assertEqual(groups, {"opu_permease": ["opuAB", "opuBB", "opuCB"]})
        self.assertEqual(set(excluded), NOT_REPORTED_V11)
        self.assertTrue(all(v == "none" for v in excluded.values()))
        self.assertFalse(set(families) & (NOT_REPORTED_V11 | {"opuAB", "opuBB", "opuCB"}))
        self.assertIn("galE", families)          # reported, though out of profile's scope

    def test_unit_both_family_and_group_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "levels.tsv"
            p.write_text("# c\nunit\treport_level\tmembers\n"
                         "x\tgene\tx\nx\tgroup\ta,b\n")
            with self.assertRaises(ValueError):
                F.read_report_levels(p)


class GroupFeatures(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.prefix = os.path.join(self.d.name, "g")
        write_gene_counts(self.prefix, {"opuAB": 1, "opuCB": 2, "ectA": 1, "mrpD": 3})

    def tearDown(self):
        self.d.cleanup()

    def test_group_is_the_sum_of_its_members(self):
        f = F.features_from_annotate(self.prefix, ["ectA", "teaA"], [],
                                     family_groups={"opu_permease": ["opuAB", "opuBB", "opuCB"]})
        self.assertEqual(f["opu_permease"], 3.0)
        self.assertNotIn("opuAB", f)
        self.assertNotIn("mrpD", f)               # not reported: not read at all
        self.assertEqual(f["teaA"], 0.0)
        # ectA, opuAB, opuCB present; mrpD is not a schema family
        self.assertEqual(f["n_families_present"], 3.0)
        self.assertEqual(f["total_osmo_copies"], 4.0)

    def test_without_groups_unchanged(self):
        a = F.features_from_annotate(self.prefix, ["ectA", "mrpD"], [])
        b = F.features_from_annotate(self.prefix, ["ectA", "mrpD"], [], family_groups={})
        self.assertEqual(a, b)
        self.assertEqual(a["total_osmo_copies"], 4.0)

    def test_feature_names_order(self):
        names = F.feature_names(["ectA"], ["ect"], ("base", "families"), {"opu_permease": ["opuAB"]})
        self.assertEqual(names, ["total_proteins", "gc", "n_contigs", "genome_mb", "ectA", "opu_permease"])


if __name__ == "__main__":
    unittest.main()
