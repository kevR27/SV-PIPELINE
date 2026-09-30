#!/usr/bin/env python3
import ast
import unittest
from pathlib import Path


class TestCandidateLocusLayout(unittest.TestCase):
    def test_plot_script_parses(self):
        path = Path(__file__).resolve().parents[1] / "plots" / "plot_candidate_locus.py"
        ast.parse(path.read_text(encoding="utf-8"))

    def test_no_tight_layout_call(self):
        path = Path(__file__).resolve().parents[1] / "plots" / "plot_candidate_locus.py"
        text = path.read_text(encoding="utf-8")
        self.assertNotIn("tight_layout(", text)
        self.assertIn("ax_methyl", text)


if __name__ == "__main__":
    unittest.main()
