#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "plots"))

from run_samplot_candidates import safe_name, short_gene_label


class TestSamplotFilenames(unittest.TestCase):
    def test_gene_label_is_short_for_large_multigene_sv(self):
        genes = ";".join(f"GENE{i}" for i in range(2000))
        label = short_gene_label("OPA1", genes)
        self.assertEqual(label, "OPA1_plus1999genes")
        self.assertLess(len(label), 50)

    def test_safe_name_is_truncated(self):
        value = "A" * 500
        name = safe_name(value, 45)
        self.assertEqual(len(name), 45)


if __name__ == "__main__":
    unittest.main()
