#!/usr/bin/env python3
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from sv_gene_effects import get_functional_context, get_sv_gene_effect


def make_row(svtype, start, end, tx_start, tx_end, location="intron", chrom="chr1"):
    return {
        "SVTYPE": svtype,
        "CHROM": chrom,
        "START": str(start),
        "END": str(end),
        "CHR2": chrom,
        "POS2": str(end),
        "ANNOTSV_GENE_ROWS_JSON": json.dumps([
            {
                "Gene_name": "TEST",
                "SV_chrom": chrom,
                "Tx_start": str(tx_start),
                "Tx_end": str(tx_end),
                "Location": location,
            }
        ]),
    }


class TestSVGeneEffects(unittest.TestCase):
    def test_gene_only_inside_inversion(self):
        effect, distance = get_sv_gene_effect(
            make_row("INV", 100, 100000, 40000, 60000),
            "TEST",
        )
        self.assertEqual(effect, "GENE_FULLY_SPANNED_BY_INVERSION")
        self.assertEqual(distance, 39900)
        self.assertEqual(
            get_functional_context("INV", effect),
            "GENE_ORIENTATION_CHANGED_COPY_NEUTRAL_REGULATORY_CONTEXT",
        )

    def test_inversion_breakpoint_in_exon(self):
        effect, distance = get_sv_gene_effect(
            make_row("INV", 450, 1000, 400, 600, location="exon"),
            "TEST",
        )
        self.assertEqual(effect, "INVERSION_BREAKPOINT_IN_EXON")
        self.assertEqual(distance, 0)
        self.assertEqual(
            get_functional_context("INV", effect),
            "DIRECT_TRANSCRIPT_DISRUPTION_POSSIBLE",
        )

    def test_inversion_breakpoint_in_intron(self):
        effect, distance = get_sv_gene_effect(
            make_row("INV", 450, 1000, 400, 600, location="intron"),
            "TEST",
        )
        self.assertEqual(effect, "INVERSION_BREAKPOINT_IN_INTRON")
        self.assertEqual(distance, 0)

    def test_inversion_breakpoint_near_gene(self):
        effect, distance = get_sv_gene_effect(
            make_row("INV", 350, 1000, 400, 600),
            "TEST",
            near_breakpoint_bp=100,
        )
        self.assertEqual(effect, "INVERSION_BREAKPOINT_NEAR_GENE")
        self.assertEqual(distance, 50)

    def test_whole_gene_deletion(self):
        effect, _ = get_sv_gene_effect(
            make_row("DEL", 100, 1000, 400, 600),
            "TEST",
        )
        self.assertEqual(effect, "WHOLE_GENE_DELETION")


if __name__ == "__main__":
    unittest.main()
