#!/usr/bin/env python3
import ast
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TestAllPythonSyntax(unittest.TestCase):
    def test_scripts_and_plots_parse(self):
        files = sorted((ROOT / "scripts").glob("*.py")) + sorted((ROOT / "plots").glob("*.py"))
        self.assertTrue(files, "No Python scripts discovered")
        for path in files:
            with self.subTest(path=str(path.relative_to(ROOT))):
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


if __name__ == "__main__":
    unittest.main()
