#!/usr/bin/env python3
"""Check the spectrum plot against the actual final candidate-table schema."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "plots"))

from build_candidate_tables import build_sv_table
import plot_gene_sv_spectrum as spectrum


class TestGeneSvSpectrum(unittest.TestCase):
    def render(self, table):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "events.tsv"
            prefix = Path(folder) / "spectrum"
            table.to_csv(source, sep="\t", index=False)
            save = spectrum.save_figure
            with patch.object(sys, "argv", [
                "plot_gene_sv_spectrum.py", "--input", str(source),
                "--out-prefix", str(prefix), "--top-n", "25",
            ]), patch.object(
                spectrum, "save_figure",
                side_effect=lambda fig, path: save(fig, path, dpi=80),
            ):
                spectrum.main()
            for extension in ["pdf", "svg", "png"]:
                self.assertGreater(prefix.with_suffix("." + extension).stat().st_size, 0)
            return pd.read_csv(str(prefix) + "_summary.tsv", sep="\t")

    def test_final_table_counts_sizes_and_breakpoints(self):
        events = pd.DataFrame([
            {
                "SV_ID": str(i), "GENES": "OPA1", "CHROM": "chr3",
                "START": "100", "END": "200", "SVTYPE": svtype,
                "SVLEN": size, "SV_EVENT_SPAN_BP": size,
                "PANEL_STATUS": "PANEL_GENE",
                "SV_GENE_RELATIONSHIP": "PARTIAL_GENE_OVERLAP",
            }
            for i, (svtype, size) in enumerate([
                ("DEL", -99999), ("DUP", 100000), ("INV", 1000000),
                ("DEL", -10000000), ("BND", "."),
            ])
        ])
        compact = build_sv_table(events, 10000)
        self.assertNotIn("SV_EVENT_SIZE_CLASS", compact)
        self.assertNotIn("BREAKPOINT_DEFINED_EVENT", compact)
        # Unknown length stays visible; repeated annotation rows count once.
        unknown = compact.iloc[[0]].copy()
        unknown["SV_ID"] = "unknown"
        unknown["SV_SIZE_GROUP"] = "UNKNOWN"
        table = pd.concat([compact, compact.iloc[[0]], unknown], ignore_index=True)
        summary = self.render(table).iloc[0]
        self.assertEqual(summary["unique_SVs"], 6)
        for category, _ in spectrum.SIZE_GROUPS:
            self.assertEqual(summary[category], 1)
        self.assertEqual(summary["breakpoint_defined_INV_BND"], 2)

    def test_legacy_table_and_empty_selection(self):
        legacy = pd.DataFrame([
            {"SV_ID": "inv", "GENE": "OPA1", "SVTYPE": "INV",
             "PANEL_STATUS": "PANEL_GENE", "SV_EVENT_SIZE_CLASS": "1MB_10MB"},
            {"SV_ID": "tra", "GENE": "TEST", "SVTYPE": "TRA",
             "PANEL_STATUS": "NONPANEL_GENE",
             "SV_EVENT_SIZE_CLASS": "BREAKEND_INTERCHROM_OR_UNRESOLVED"},
        ])
        summary = self.render(legacy).set_index("gene")
        self.assertEqual(summary.loc["OPA1", "1MB_10MB"], 1)
        self.assertEqual(summary.loc["TEST", "BREAKEND"], 1)
        self.assertEqual(summary["breakpoint_defined_INV_BND"].sum(), 2)
        legacy["GENE"] = "."
        self.assertTrue(self.render(legacy).empty)


if __name__ == "__main__":
    unittest.main()
