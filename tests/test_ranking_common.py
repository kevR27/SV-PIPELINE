#!/usr/bin/env python3
import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ranking_common import (
    gene_relevance,
    inferred_effect_class,
    mechanism_inheritance_summary,
    semantic_engine,
)


class TestRankingCommon(unittest.TestCase):
    def test_asymmetric_query_coverage_does_not_penalize_extra_gene_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            edges = Path(tmp) / "edges.tsv"
            rows = [
                {
                    "subject": "HGNC:A",
                    "object": "HP:0000001",
                    "predicate": "biolink:has_phenotype",
                    "category": "biolink:GeneToPhenotypicFeatureAssociation",
                },
                {
                    "subject": "HGNC:A",
                    "object": "HP:0009999",
                    "predicate": "biolink:has_phenotype",
                    "category": "biolink:GeneToPhenotypicFeatureAssociation",
                },
                {
                    "subject": "HGNC:B",
                    "object": "HP:0000002",
                    "predicate": "biolink:has_phenotype",
                    "category": "biolink:GeneToPhenotypicFeatureAssociation",
                },
                {
                    "subject": "HGNC:C",
                    "object": "HP:0000003",
                    "predicate": "biolink:has_phenotype",
                    "category": "biolink:GeneToPhenotypicFeatureAssociation",
                },
            ]
            with edges.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(
                    fh,
                    fieldnames=["subject", "object", "predicate", "category"],
                    delimiter="\t",
                )
                writer.writeheader()
                writer.writerows(rows)

            engine = semantic_engine(str(edges))
            query = {"HP:0000001"}
            pleiotropic_gene = {"HP:0000001", "HP:0009999"}

            self.assertAlmostEqual(
                engine["normalized_query_coverage"](query, pleiotropic_gene),
                1.0,
                places=6,
            )
            self.assertLess(
                engine["normalized_bma"](pleiotropic_gene, query),
                1.0,
            )

    def test_partial_duplication_is_possible_disruption_not_pure_copy_gain(self):
        self.assertEqual(
            inferred_effect_class("DUP", "PARTIAL_GENE_OVERLAP"),
            "DISRUPTION_OR_LOF",
        )
        self.assertEqual(
            inferred_effect_class("DUP", "BREAKPOINT_WITHIN_TRANSCRIPT"),
            "DISRUPTION_OR_LOF",
        )
        self.assertEqual(
            inferred_effect_class("DUP", "WHOLE_GENE_DOSAGE_CONTEXT"),
            "COPY_GAIN",
        )

    def test_partial_duplication_can_match_dominant_lof_model(self):
        result = mechanism_inheritance_summary(
            {
                "GENE": "TEST",
                "CHROM": "chr1",
                "SVTYPE": "DUP",
                "SV_GENE_RELATIONSHIP": "PARTIAL_GENE_OVERLAP",
                "GENCC_MOI": "Autosomal dominant",
                "OMIM_INHERITANCE": ".",
                "ALLELE_GENOTYPE": "HET",
                "CLINGEN_HI": ".",
                "CLINGEN_TS": ".",
            }
        )
        self.assertEqual(result["effect_class"], "DISRUPTION_OR_LOF")
        self.assertEqual(
            result["category"],
            "SUPPORTED_DISEASE_MECHANISM",
        )

    def test_gene_relevance_tiers_are_reachable_on_normalized_scale(self):
        self.assertEqual(gene_relevance(0.80, 2.0)["tier"], "HIGH")
        self.assertEqual(gene_relevance(0.50, 0.0)["tier"], "MODERATE")
        self.assertEqual(gene_relevance(0.25, 0.0)["tier"], "SUPPORTING")
        self.assertEqual(gene_relevance(0.10, 0.0)["tier"], "LIMITED")


if __name__ == "__main__":
    unittest.main()
