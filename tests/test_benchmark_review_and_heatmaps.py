"""Regression tests for per-patient benchmark review and readable heatmaps."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "plots"))
from summarize_benchmarking_review import report
from plot_gene_hpo_heatmap import make_heatmap
from plot_candidate_evidence_matrix import render_matrix
from candidate_review import domains_for, STATE_COLORS


class BenchmarkOverviewTests(unittest.TestCase):
    def test_stability_is_not_reported_as_truth_benchmark(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            compare = folder / "comparison.json"
            compare.write_text(json.dumps({
                "sample_id": "11446", "top_k": {
                    "5": {"custom_n": 5, "genomiser_n": 5, "jaccard": 0.25}
                }
            }))
            stability = folder / "stability.tsv"
            stability.write_text(
                "scenario\tstatus\tbaseline_top_n_retained\n"
                "baseline\tEVALUATED\t1.0\n"
                "at_least_2_callers\tEVALUATED\t0.5\n"
            )
            rows, manifest = report("11446", compare, stability=stability)
            cats = {r["CATEGORY"] for r in rows}
            self.assertEqual(cats, {
                "METHOD_AGREEMENT", "SENSITIVITY_ANALYSIS",
                "CONFIRMED_TRUTH_BENCHMARK",
            })
            self.assertEqual(manifest["confirmed_truth_benchmark"], "NOT_CONFIGURED")
            self.assertEqual(manifest["ranking_stability"], "EVALUATED")
            self.assertEqual(
                next(x for x in rows if x["CATEGORY"] == "CONFIRMED_TRUTH_BENCHMARK")["STATUS"],
                "NOT_CONFIGURED",
            )
            self.assertEqual(
                next(x for x in rows if x["SCENARIO"] == "at_least_2_callers")["VALUE"],
                "0.5",
            )

    def test_truth_file_is_checked_against_sample(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            comparison = folder / "comparison.json"
            comparison.write_text(json.dumps({"sample_id": "11446", "top_k": {}}))
            truth = folder / "truth.json"
            truth.write_text(json.dumps({"sample_id": "other", "ground_truth_records": 1}))
            with self.assertRaisesRegex(ValueError, "another sample"):
                report("11446", comparison, truth=truth)
            truth.write_text(json.dumps({"sample_id": "11446", "ground_truth_records": 2}))
            rows, manifest = report("11446", comparison, truth=truth)
            self.assertEqual(manifest["confirmed_truth_benchmark"], "EVALUATED")
            self.assertEqual(rows[-1]["VALUE"], 2)


class HeatmapLayoutTests(unittest.TestCase):
    def test_hpo_heatmap_expands_and_rotates_labels(self):
        mat = pd.DataFrame(np.zeros((10, 14)), index=[f"G{i}" for i in range(10)])
        labels = [f"HP:000{i:04d}\nLong phenotype\nlabel number {i}" for i in range(14)]
        state = {}
        def capture(fig, *args, **kwargs):
            state["width"], state["height"] = fig.get_size_inches()
            state["angles"] = [lab.get_rotation() for lab in fig.axes[0].get_xticklabels()]
            return []
        with patch("plot_gene_hpo_heatmap.save_figure", side_effect=capture):
            make_heatmap(mat, labels, set(), Path("/tmp/unused"), "HPO", list(range(14)))
        self.assertGreater(state["width"], 15)
        self.assertGreater(state["height"], 7)
        self.assertEqual(set(state["angles"]), {90.0})

    def test_candidate_matrix_allocates_room_per_row_and_domain(self):
        domains = domains_for("lrs")
        first_state = next(iter(STATE_COLORS))
        observed = {d: (first_state, "Evidence\nunder review") for d in domains}
        data = pd.DataFrame([
            {"SV_ID": f"sv{i}", "GENE": f"G{i}", "CHROM": "chr1",
             "START": str(i * 1000), "SVTYPE": "DEL", "_rank": float(i + 1)}
            for i in range(12)
        ])
        metrics = {}
        def capture(fig, *args, **kwargs):
            metrics["width"], metrics["height"] = fig.get_size_inches()
            metrics["rows"] = len(fig.axes[0].get_yticklabels())
            return []
        with patch("plot_candidate_evidence_matrix.event_evidence", return_value=observed):
            with patch("plot_candidate_evidence_matrix.save_figure", side_effect=capture):
                render_matrix(data, Path("/tmp/unused"), "Evidence", platform="lrs")
        self.assertGreaterEqual(metrics["width"], len(domains) * 2.15)
        self.assertGreaterEqual(metrics["height"], 12 * 0.88)
        self.assertEqual(metrics["rows"], 12)


if __name__ == "__main__":
    unittest.main()
