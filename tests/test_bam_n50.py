#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from calc_bam_n50 import n50


class TestBamN50(unittest.TestCase):
    def test_n50_is_base_weighted(self):
        # Total bases = 100 + 80 + 20 = 200; half = 100, so N50 = 100.
        self.assertEqual(n50([100, 80, 20]), 100)

    def test_n50_not_median(self):
        # Median is 10, but the 100 bp read alone contains >50% of bases.
        self.assertEqual(n50([100, 10, 10, 10, 10]), 100)

    def test_empty(self):
        self.assertEqual(n50([]), 0)


if __name__ == "__main__":
    unittest.main()
