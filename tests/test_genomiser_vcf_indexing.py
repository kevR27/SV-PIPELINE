"""Regression tests: Genomiser must index its local VCFs before concat."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_genomiser import prepare_vcf


class GenomiserIndexTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.folder = Path(tmp.name)
        self.source = self.folder / "original.vcf.gz"
        self.source.write_bytes(b"placeholder")
        self.target = self.folder / "genomiser" / "small.renamed.vcf.gz"
        self.target.parent.mkdir()
        self.namefile = self.folder / "genomiser" / "sample_name.txt"

    def test_reheader_indexes_new_vcf_before_reading(self):
        with patch("run_genomiser.samples", side_effect=[["old_id"], ["19545"]]), \
             patch("run_genomiser.run") as execute:
            prepared, old = prepare_vcf(self.source, "19545", self.target, self.namefile)
        self.assertEqual(prepared, str(self.target))
        self.assertEqual(old, "old_id")
        self.assertEqual(self.namefile.read_text(), "19545\n")
        self.assertEqual(execute.call_args_list, [
            call(["bcftools", "reheader", "-s", self.namefile,
                  "-o", self.target, self.source]),
            call(["bcftools", "index", "-f", "-t", str(self.target)]),
        ])

    def test_matching_sample_and_existing_index_reuses_input(self):
        Path(str(self.source) + ".tbi").write_bytes(b"index")
        with patch("run_genomiser.samples", return_value=["19545"]), \
             patch("run_genomiser.run") as execute:
            prepared, old = prepare_vcf(self.source, "19545", self.target, self.namefile)
        self.assertEqual(prepared, str(self.source))
        self.assertEqual(old, "19545")
        execute.assert_not_called()

    def test_missing_original_index_uses_local_symlink(self):
        with patch("run_genomiser.samples", side_effect=[["19545"], ["19545"]]), \
             patch("run_genomiser.run") as execute:
            prepared, old = prepare_vcf(self.source, "19545", self.target, self.namefile)
        self.assertEqual(prepared, str(self.target))
        self.assertEqual(old, "19545")
        self.assertTrue(self.target.is_symlink())
        self.assertEqual(self.target.resolve(), self.source.resolve())
        execute.assert_called_once_with(
            ["bcftools", "index", "-f", "-t", str(self.target)]
        )


if __name__ == "__main__":
    unittest.main()
