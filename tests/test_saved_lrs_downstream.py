"""Tests for archived LRS source discovery and BAM-free reanalysis contracts."""
import tempfile
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_saved_lrs_outputs import discover_samples, required_sources, optional_sources, audit


class SavedLrsDownstreamTests(unittest.TestCase):
    def test_discovery_uses_sample_folder_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ["21930", "11446", "cohort_analysis", "sample_A", "25428"]:
                (root / name).mkdir()
            self.assertEqual(discover_samples(root), ["11446", "21930", "25428"])

    def test_sample_11446_resolves_original_relative_paths(self):
        required = required_sources("11446")
        self.assertEqual(
            required["merged_structural_variants"],
            "11446/sv/merged/11446_merged_SV.vcf.gz",
        )
        self.assertEqual(
            required["small_variants"], "11446/snp_clair3/11446.vcf.gz"
        )
        self.assertIn("11446/coverage/binned/11446.10kb.regions.bed.gz",
                      required.values())
        self.assertFalse(any(".bam" in item for item in required.values()))

    def test_audit_classifies_saved_and_missing_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            small = root / required_sources("11446")["small_variants"]
            small.parent.mkdir(parents=True)
            small.write_bytes(b"dummy")
            rows = audit(root, "11446")
            self.assertEqual(
                [row["status"] for row in rows if row["component"] == "small_variants"],
                ["OK"],
            )
            self.assertTrue(any(row["category"] == "REQUIRED" and
                                row["status"] == "MISSING" for row in rows))
            self.assertEqual(len(rows), len(required_sources("11446")) +
                             len(optional_sources("11446")))

    def test_bam_free_postprocess_hides_read_dependent_rules(self):
        original = (ROOT / "snakemake_pipelines" / "lrs" /
                    "Snakefile_LRS_postprocess").read_text(encoding="utf-8")
        extensions = (ROOT / "snakemake_pipelines" / "lrs" /
                      "thesis_extensions.smk").read_text(encoding="utf-8")
        snakefile = (ROOT / "snakemake_pipelines" /
                     "lrs_downstream_reanalysis" / "Snakefile").read_text(encoding="utf-8")
        self.assertIn("bam_free_downstream", snakefile)
        self.assertNotIn("samplot_candidate_breakpoints", original)
        self.assertIn("if not BAM_FREE_DOWNSTREAM:", extensions)
        self.assertIn('coverage=PATH + "{sample}/coverage/binned/', extensions)
        self.assertIn('include: "../lrs/Snakefile_LRS_postprocess"', snakefile)
        self.assertNotIn("rule binned_read_depth:", snakefile)
        self.assertNotIn("rule samplot_candidate_breakpoints:", snakefile)


if __name__ == "__main__":
    unittest.main()
