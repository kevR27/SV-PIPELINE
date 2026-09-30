#!/usr/bin/env python3
import ast
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TestThesisExtensionScripts(unittest.TestCase):
    def test_new_python_scripts_parse(self):
        for relative in [
            "scripts/build_hon_hpo_reference.py",
            "scripts/hpo_semantic_similarity.py",
            "scripts/integrate_depth_into_candidates.py",
            "scripts/build_large_sv_gene_context.py",
            "scripts/build_population_gene_effect_table.py",
        ]:
            with self.subTest(relative=relative):
                path = ROOT / relative
                ast.parse(path.read_text(encoding="utf-8"))

    def test_hpo_semantic_similarity_has_safe_no_patient_mode(self):
        text = (ROOT / "scripts/hpo_semantic_similarity.py").read_text(encoding="utf-8")
        self.assertIn("PATIENT_HPO_NOT_AVAILABLE", text)
        self.assertIn("gene-phenotypes", text)

    def test_large_sv_gene_context_reads_genes_from_bed_interval(self):
        text = (ROOT / "scripts/build_large_sv_gene_context.py").read_text(encoding="utf-8")
        self.assertIn('bed["GENE_END"].gt(start)', text)
        self.assertIn('bed["GENE_START"].lt(end)', text)


if __name__ == "__main__":
    unittest.main()
