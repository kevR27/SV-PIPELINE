"""Regression coverage for Genomiser-only large deletion representation."""
import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_genomiser_sv_vcf import prepare_sv_vcf


class SymbolicDeletionTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.source = self.root / "original.sv.vcf.gz"
        self.output = self.root / "genomiser.sv.vcf.gz"
        self.audit = self.root / "symbolic_deletions.tsv"

    def write_source(self, records):
        with gzip.open(self.source, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n')
            fh.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">\n')
            fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t19545\n")
            for record in records:
                fh.write(record + "\n")

    def fake_bcftools(self, command):
        if command[:2] == ["bcftools", "view"]:
            source = Path(command[-1])
            with source.open("rt") as fh:
                contents = fh.read()
            with gzip.open(self.output, "wt") as fh:
                fh.write(contents)
        elif command[:2] == ["bcftools", "index"]:
            Path(str(self.output) + ".tbi").write_text("test-index")
        else:
            raise AssertionError(f"Unexpected command {command}")

    def test_convert_large_deletion_but_not_bnd_or_small_indel(self):
        long_ref = "A" + "C" * 1200
        self.write_source([
            f"chr1\t100\tbig_del\t{long_ref}\tA\t60\tPASS\tSVTYPE=DEL;END=1300;SVLEN=-1200\tGT:GQ\t0/1:60",
            "chr1\t180\tbnd\tC\tC]chr5:25]\t60\tPASS\tSVTYPE=BND;END=25\tGT\t0/1",
            "chr1\t250\tsmall\tAC\tA\t60\tPASS\t.\tGT\t0/1",
        ])
        with patch("prepare_genomiser_sv_vcf._run", side_effect=self.fake_bcftools):
            summary = prepare_sv_vcf(self.source, self.output, self.audit)
        with gzip.open(self.output, "rt") as fh:
            records = [line.rstrip("\n") for line in fh if not line.startswith("#")]
            header = []
        first = records[0].split("\t")
        self.assertEqual(first[3:5], ["A", "<DEL>"])
        self.assertIn("END=1300", first[7])
        self.assertIn("SVLEN=-1200", first[7])
        self.assertEqual(first[9], "0/1:60")
        self.assertEqual(records[1].split("\t")[7], "SVTYPE=BND;END=25")
        self.assertEqual(records[2].split("\t")[3:5], ["AC", "A"])
        self.assertEqual(summary["symbolic_deletions"], 1)
        with self.audit.open(newline="") as handle:
            audit = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(audit), 1)
        self.assertEqual(audit[0]["ID"], "big_del")
        self.assertEqual(audit[0]["OLD_REF_LENGTH"], "1201")
        with gzip.open(self.source, "rt") as handle:
            self.assertIn(long_ref, handle.read())
        self.assertFalse(self.output.with_suffix("").exists())

    def test_detect_inconsistent_old_interval_in_audit(self):
        long_ref = "G" + "T" * 1200
        self.write_source([
            f"chr2\t1000\tconflict\t{long_ref}\tG\t60\tPASS\tSVTYPE=DEL;END=1001;SVLEN=-500\tGT\t0/1"
        ])
        with patch("prepare_genomiser_sv_vcf._run", side_effect=self.fake_bcftools):
            prepare_sv_vcf(self.source, self.output, self.audit)
        with self.audit.open(newline="") as handle:
            row = list(csv.DictReader(handle, delimiter="\t"))[0]
        self.assertEqual(row["END_LENGTH_DISAGREEMENT"], "True")
        self.assertEqual(row["SVLEN_DISAGREEMENT"], "True")
        self.assertEqual(row["NEW_END"], "2200")

    def test_reject_unconvertible_megabase_sequence(self):
        self.write_source([
            f"chr4\t20\tcomplex\t{'A' * 1100}\tC\t40\tPASS\tSVTYPE=DEL\tGT\t0/1"
        ])
        with patch("prepare_genomiser_sv_vcf._run", side_effect=self.fake_bcftools):
            with self.assertRaisesRegex(ValueError, "Unconverted explicit allele"):
                prepare_sv_vcf(self.source, self.output, self.audit, maximum_unsupported_length=1000)


if __name__ == "__main__":
    unittest.main()
