"""Unit tests for catalog_benchmark.py's per-read rebuild of profile --catalog.

catalog_benchmark.py and compute_metrics.py need pandas (and compute_metrics.py
numpy, matplotlib, scikit-learn, pyyaml: the metrics container). Run from the
repository root:  python -m unittest discover -s tests/python -v
"""
import importlib.util
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


try:
    import pandas as pd
    cm = load("compute_metrics", REPO / "modules/local/benchmark/metrics/templates/compute_metrics.py")
    cb = load("catalog_benchmark", REPO / "modules/local/benchmark/catalogscore/templates/catalog_benchmark.py")
except ImportError:
    cm = cb = None


@unittest.skipIf(cb is None, "catalog_benchmark.py's packages are not installed")
class HybridCalls(unittest.TestCase):
    """--catalog_only settles an ambiguous read; profile with translated search
    (osmotool >= 0.12.0) sends it on, so the rebuild gives it its translated call."""
    def setUp(self):
        self.settled = pd.DataFrame([
            {"read_id": "r1/1", "sample": "s", "decision": "call", "label": "proP"},
            {"read_id": "r2/1", "sample": "s", "decision": "drop", "label": ""},
            {"read_id": "r3/1", "sample": "s", "decision": "ambiguous", "label": "kdpA,drop"},
        ])
        self.translated = pd.DataFrame([
            {"read_id": rid, "sample": "s", "family": fam, "bitscore": 50.0, "evalue": 1e-5}
            for rid, fam in [("r1/1", "proX"), ("r2/1", "proP"), ("r3/1", "kdpA"), ("r4/1", "nhaA")]
        ])

    def calls(self, fallback):
        out = cb.hybrid_calls(cm, self.translated, self.settled, fallback=fallback)
        return dict(zip(out["read_id"], out["family"]))

    def test_with_fallback(self):
        # r1 the catalog's call, r2 dropped, r3 ambiguous -> translated, r4 unsettled -> translated
        self.assertEqual(self.calls(True), {"r1/1": "proP", "r3/1": "kdpA", "r4/1": "nhaA"})

    def test_catalog_only(self):
        # no translated search: an ambiguous read is no call
        self.assertEqual(self.calls(False), {"r1/1": "proP"})


if __name__ == "__main__":
    unittest.main()
