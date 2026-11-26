"""Unit tests for real_read_truth.py's mate placement (standard library only).

Run from the repository root:  python -m unittest discover -s tests/python -v
"""
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "modules/local/validation/readrecall/templates/real_genome_validation"))
import real_read_truth as rt  # noqa: E402


def revcomp(s: str) -> str:
    return s.translate(rt.COMPLEMENT)[::-1]


class ReadStartTest(unittest.TestCase):
    """wgsim names carry the fragment's outer coordinates; either mate can be at either end."""
    def setUp(self):
        rng = random.Random(1)
        self.genome = "".join(rng.choice("ACGT") for _ in range(2000))
        self.frag_start, self.frag_end, self.n = 301, 800, 150   # 1-based, inclusive

    def left(self):
        return self.genome[self.frag_start - 1:self.frag_start - 1 + self.n]

    def right(self):
        return revcomp(self.genome[self.frag_end - self.n:self.frag_end])

    def test_left_end_read_forward(self):
        self.assertEqual(rt.read_start(self.left(), self.genome, self.frag_start, self.frag_end), 301)

    def test_right_end_read_reverse_complement(self):
        self.assertEqual(rt.read_start(self.right(), self.genome, self.frag_start, self.frag_end), 651)

    def test_sequencing_errors_do_not_move_the_read(self):
        read = list(self.right())
        for i in (3, 40, 77, 120, 149):
            read[i] = "A" if read[i] != "A" else "C"
        self.assertEqual(rt.read_start("".join(read), self.genome, self.frag_start, self.frag_end), 651)

    def test_mate_number_is_not_used(self):
        # the old rule put mate 2 at frag_end (after the fragment); placement depends on sequence only
        m = rt.HEADER_RE.match("chr_301_800_1:0:0_2:0:0_1a/2")
        self.assertEqual((m.group("frag_start"), m.group("frag_end")), ("301", "800"))


class OverlapTest(unittest.TestCase):
    def test_min_overlap(self):
        regions = [("proP", "chr", 700, 1500)]
        self.assertEqual(rt.overlapping_family("chr", 651, 800, regions, 20), "proP")
        self.assertIsNone(rt.overlapping_family("chr", 301, 450, regions, 20))
        self.assertIsNone(rt.overlapping_family("chr", 532, 681, regions, 20))   # 19 bp overlap


if __name__ == "__main__":
    unittest.main()
