#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_population_gene_effect_table import (
    benign_afmax_class,
    gnomad_context,
    needlr_class,
    population_class,
)


class TestPopulationGeneEffect(unittest.TestCase):
    def test_needlr_frequency_classes(self):
        self.assertEqual(
            needlr_class({"NEEDLR_AF": "0", "POPULATION_STATUS": "RARE"}),
            "NOT_OBSERVED_IN_NEEDLR_CONTROLS",
        )
        self.assertEqual(
            needlr_class({"NEEDLR_AF": "0.005", "POPULATION_STATUS": "RARE"}),
            "RARE_NEEDLR_LE_0.01",
        )
        self.assertEqual(
            needlr_class({"NEEDLR_AF": "0.2", "POPULATION_STATUS": "COMMON"}),
            "VERY_COMMON_NEEDLR_GE_0.05",
        )

    def test_no_needlr_match_is_not_af_zero(self):
        self.assertEqual(
            needlr_class({"NEEDLR_AF": ".", "POPULATION_STATUS": "NO_POPULATION_MATCH"}),
            "NO_NEEDLR_MATCH_AF_UNKNOWN",
        )

    def test_annotsv_afmax_is_benign_region_context(self):
        row = {
            "NEEDLR_AF": ".",
            "POPULATION_STATUS": "UNKNOWN",
            "ANNOTSV_BENIGN_AFMAX": "0.15",
            "GNOMAD_SV_OVERLAP": "YES",
            "ANNOTSV_BENIGN_DB_SOURCE": "gnomAD;DGV",
        }
        self.assertEqual(
            benign_afmax_class(row),
            "BENIGN_REGION_AFMAX_GE_0.05",
        )
        self.assertEqual(
            gnomad_context(row),
            "GNOMAD_INCLUDED_IN_ANNOTSV_BENIGN_OVERLAP",
        )
        self.assertEqual(
            population_class(row),
            "COMMON_BENIGN_REGION_OVERLAP_CONTEXT",
        )

    def test_exact_gnomad_frequency_contributes(self):
        row = {
            "NEEDLR_AF": ".",
            "POPULATION_STATUS": "UNKNOWN",
            "GNOMAD_SV_EXACT_MATCH": "YES",
            "GNOMAD_SV_AF": "0.0005",
            "ANNOTSV_BENIGN_AFMAX": ".",
        }
        self.assertEqual(
            population_class(row),
            "LOW_FREQUENCY_BY_EXACT_GNOMAD",
        )

    def test_frequency_conflict_is_preserved(self):
        row = {
            "NEEDLR_AF": "0.2",
            "POPULATION_STATUS": "COMMON",
            "GNOMAD_SV_EXACT_MATCH": "YES",
            "GNOMAD_SV_AF": "0.0001",
            "ANNOTSV_BENIGN_AFMAX": ".",
        }
        self.assertEqual(
            population_class(row),
            "FREQUENCY_SOURCES_CONFLICT",
        )

    def test_low_frequency_is_not_named_pathogenic(self):
        cls = population_class(
            {
                "NEEDLR_AF": "0.0001",
                "POPULATION_STATUS": "RARE",
                "ANNOTSV_BENIGN_AFMAX": ".",
            }
        )
        self.assertEqual(cls, "LOW_FREQUENCY_BY_NEEDLR")
        self.assertNotIn("PATHOGEN", cls)
        self.assertNotIn("BENIGN", cls)


if __name__ == "__main__":
    unittest.main()
