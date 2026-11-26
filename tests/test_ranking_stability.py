"""Small regression checks for the descriptive SV-gene shortlist sensitivity."""
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_ranking_stability import analyze


class RankingStabilityTests(unittest.TestCase):
    def make_table(self, folder, data):
        path = Path(folder) / "ranked.tsv"
        pd.DataFrame(data).to_csv(path, sep="\t", index=False)
        return path

    def test_population_and_caller_scenarios_do_not_rescore_events(self):
        with tempfile.TemporaryDirectory() as folder:
            ranked = self.make_table(folder, {
                "SV_ID": ["s1", "s2", "s3"],
                "GENE": ["OPA1", "DNAJC30", "MFN2"],
                "FINAL_EVENT_RANK_GLOBAL": [1, 2, 3],
                "CALLER_COUNT": [1, 2, 2],
                "GNOMAD_SV_AF": [0.001, 0.2, 0.8],
                "GNOMAD_SV_EXACT_MATCH": ["YES", "YES", "NO"],
            })
            summary, details, n_events, n_genes = analyze(ranked, top_n=2)
            result = summary.set_index("scenario")
            self.assertEqual((n_events, n_genes), (3, 3))
            self.assertEqual(result.loc["baseline", "baseline_top_n_retained"], 1.0)
            # MFN2's unmatched population AF must remain unknown, not common.
            self.assertEqual(result.loc["af_le_0.01_unknown_review", "shared_top_n"], 1)
            self.assertEqual(result.loc["af_le_0.01_unknown_review", "unknown_af_events_retained"], 1)
            self.assertEqual(result.loc["at_least_2_callers", "shared_top_n"], 1)
            self.assertEqual(
                details.loc[(details.scenario == "at_least_2_callers") &
                            (details.gene == "DNAJC30"), "scenario_gene_rank"].iloc[0], 1
            )

    def test_missing_caller_counts_are_not_assumed_zero(self):
        with tempfile.TemporaryDirectory() as folder:
            ranked = self.make_table(folder, {
                "SV_ID": ["s1", "s2"],
                "GENE": ["OPA1", "MFN2"],
                "FINAL_EVENT_RANK_GLOBAL": [1, 2],
            })
            summary, _, _, _ = analyze(ranked, top_n=2)
            checks = summary.set_index("scenario")
            self.assertEqual(checks.loc["at_least_2_callers", "status"],
                             "NOT_EVALUABLE_NO_CALLER_COUNTS")
            self.assertEqual(checks.loc["baseline", "status"], "EVALUATED")


if __name__ == "__main__":
    unittest.main()
