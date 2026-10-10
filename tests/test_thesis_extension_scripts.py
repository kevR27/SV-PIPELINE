#!/usr/bin/env python3
import ast
import sys
from pathlib import Path
import tempfile
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import binom
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import calculate_xci_skew as xci
import unittest


ROOT = Path(__file__).resolve().parents[1]


class TestThesisExtensionScripts(unittest.TestCase):
    def test_new_python_scripts_parse(self):
        for relative in [
            "scripts/build_hon_hpo_reference.py",
            "scripts/hpo_semantic_similarity.py",
            "scripts/integrate_depth_into_candidates.py",
            "scripts/build_large_sv_gene_context.py",
            "scripts/build_population_gene_effect_table.py",
            "plots/plot_mitocarta_sv_genes.py",
            "plots/run_samplot_candidates.py",
            "plots/plot_candidate_locus.py",
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


    def test_core_allele_assessment_does_not_claim_exon_breakpoint(self):
        text = (ROOT / "scripts" / "assess_sv_alleles.py").read_text(encoding="utf-8")
        self.assertNotIn("EXONIC_BREAKPOINT_POSSIBLE", text)
        self.assertNotIn("INTRONIC_BREAKPOINT_POSSIBLE", text)
        self.assertIn("TRANSCRIPT_BREAKPOINT_POSSIBLE", text)


class XciNumericalTests(unittest.TestCase):
    def test_empty_block_input_reports_no_estimate(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            raw = pd.DataFrame(columns=["PS", "BLOCK_START", "BLOCK_END", "H1_Xa", "H1_Xi", "H2_Xa", "H2_Xi"])
            raw.to_csv(folder / "raw.tsv", sep="\t", index=False)
            pd.DataFrame(columns=["WHATSHAP_PS", "SHARED_PHASED_SNVS"]).to_csv(folder / "phase.tsv", sep="\t", index=False)
            pd.DataFrame([{"haplotagged_reads": 0}]).to_csv(folder / "hap.tsv", sep="\t", index=False)
            pd.DataFrame([{"WHATSHAP_CHRX_HET_SNVS": 100}]).to_csv(folder / "phase_summary.tsv", sep="\t", index=False)
            pd.DataFrame([{"chrom": "chr1", "mean": 30}, {"chrom": "chrX", "mean": 30}]).to_csv(folder / "coverage.tsv", sep="\t", index=False)
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/calculate_xci_skew.py"),
                "--block-skew", str(folder / "raw.tsv"), "--phase-blocks", str(folder / "phase.tsv"),
                "--phase-summary", str(folder / "phase_summary.tsv"), "--haplotag-summary", str(folder / "hap.tsv"),
                "--mosdepth-summary", str(folder / "coverage.tsv"), "--blocks-output", str(folder / "blocks.tsv"),
                "--summary-output", str(folder / "summary.tsv"), "--sensitivity-output", str(folder / "sensitivity.tsv"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            summary = pd.read_csv(folder / "summary.tsv", sep="\t", dtype=str)
            self.assertEqual(summary["XCI_ANALYSIS_STATUS"].iloc[0], "NO_INFORMATIVE_PHASED_XCI_BLOCKS")
            self.assertEqual(summary["GLOBAL_FOLDED_SKEW_P"].iloc[0], ".")

    def test_folded_probabilities_normalize_and_match_binomial(self):
        for n in [1, 2, 5, 30]:
            for p in [0, 0.2, 0.5]:
                probabilities = [np.exp(xci.folded_logpmf(x, n, p)) for x in range(n // 2 + 1)]
                self.assertAlmostEqual(sum(probabilities), 1)
                for x, observed in enumerate(probabilities):
                    expected = binom.pmf(x, n, p)
                    if 2*x != n:
                        expected += binom.pmf(n-x, n, p)
                    self.assertAlmostEqual(observed, expected)

    def test_exact_boundaries_and_uninformative_blocks(self):
        for x, expected in [(0, 0), (50, 0.5)]:
            block = pd.DataFrame({"SUCCESS_FOLDED": [x], "TRIALS": [100]})
            estimate = xci.fit_folded_mle(block)
            self.assertEqual(estimate, expected)
            low, high = xci.profile_likelihood_ci(block, estimate)
            self.assertTrue(low <= estimate <= high)
        one_read = pd.DataFrame({"SUCCESS_FOLDED": [0, 0], "TRIALS": [1, 1]})
        self.assertTrue(np.isnan(xci.fit_folded_mle(one_read)))

    def test_vector_likelihood_agrees_with_individual_blocks_and_grid(self):
        blocks = pd.DataFrame({"SUCCESS_FOLDED": [5, 20, 31], "TRIALS": [20, 50, 80]})
        xs, ns = blocks.SUCCESS_FOLDED.values, blocks.TRIALS.values
        for p in [0.01, 0.2, 0.5]:
            expected = -sum(xci.folded_logpmf(x, n, p) for x, n in zip(xs, ns))
            self.assertAlmostEqual(xci.negative_log_likelihood(p, xs, ns), expected)
        estimate = xci.fit_folded_mle(blocks)
        grid = np.linspace(0.001, 0.5, 1000)
        best_grid = grid[np.argmin([xci.negative_log_likelihood(p, xs, ns) for p in grid])]
        self.assertAlmostEqual(estimate, best_grid, delta=0.001)

    def test_balanced_lrt_and_threshold_equality(self):
        balanced = pd.DataFrame({"SUCCESS_FOLDED": [50], "TRIALS": [100]})
        self.assertEqual(xci.balanced_lrt(balanced, 0.5), (0, 1))
        self.assertEqual(xci.conventional_context(0.2, 0.3, 0.2, 0.1), "AT_OR_BEYOND_80_20_BUT_BELOW_90_10")


if __name__ == "__main__":
    unittest.main()
