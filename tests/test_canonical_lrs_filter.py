#!/usr/bin/env python3
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TestCanonicalLRSFilter(unittest.TestCase):
    def test_all_three_lrs_filters_request_canonical_only(self):
        snakefile = (
            ROOT / "snakemake_pipelines" / "lrs" / "Snakefile_LRS_update"
        ).read_text(encoding="utf-8")
        self.assertIn('EVIDENCE_CANONICAL_FLAG = "--canonical-only"', snakefile)
        self.assertGreaterEqual(snakefile.count("{EVIDENCE_CANONICAL_FLAG}"), 3)

    def test_canonical_set_contains_only_primary_chromosomes(self):
        script = (ROOT / "scripts" / "filter_sv_evidence.py").read_text(encoding="utf-8")
        self.assertIn('"chrX"', script)
        self.assertIn('"chrY"', script)
        self.assertIn('"chrM"', script)
        self.assertIn('range(1, 23)', script)
        self.assertIn("NON_CANONICAL_BREAKEND_CHROM", script)


if __name__ == "__main__":
    unittest.main()
