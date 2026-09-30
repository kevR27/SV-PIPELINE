#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_candidate_tables import build_gene_table


class TestFinalCandidateTables(unittest.TestCase):
    def test_gene_relevance_does_not_use_sv_count(self):
        ranking = pd.DataFrame([
            {
                "gene": "OPA1",
                "panel_gene": "YES",
                "human_HPO_count": "12",
                "optic_neuropathy_anchor_HPO_count": "2",
                "phenotype_score": "6.0",
                "gene_disease_evidence_score": "4.0",
                "gene_disease_evidence_level": "DEFINITIVE",
                "SV_count": "99",
                "SV_types": "DEL;INS",
                "SV_count_lt100kb": "50",
                "SV_count_100kb_to_1Mb": "20",
                "SV_count_1Mb_to_10Mb": "20",
                "SV_count_ge10Mb": "9",
                "breakpoint_defined_INV_BND_count": "0",
                "candidate_group": "PANEL_GENE",
            }
        ])
        out = build_gene_table(ranking, pd.DataFrame())
        self.assertEqual(float(out.iloc[0]["GENE_RELEVANCE_SCORE"]), 10.0)
        self.assertEqual(out.iloc[0]["SV_COUNT"], "99")


if __name__ == "__main__":
    unittest.main()
