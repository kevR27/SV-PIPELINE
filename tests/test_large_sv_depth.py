#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from summarize_large_sv_depth import choose_plot_bin


class TestLargeSVDepthBins(unittest.TestCase):
    def test_adaptive_plot_bins(self):
        self.assertEqual(choose_plot_bin(500_000), 50_000)
        self.assertEqual(choose_plot_bin(2_000_000), 100_000)
        self.assertEqual(choose_plot_bin(20_000_000), 250_000)
        self.assertEqual(choose_plot_bin(80_000_000), 500_000)


if __name__ == "__main__":
    unittest.main()
