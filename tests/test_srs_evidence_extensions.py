from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRS = ROOT / "scripts" / "srs"


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SRSEvidenceExtensionTests(unittest.TestCase):
    def run_script(self, script, *args):
        subprocess.run(
            [sys.executable, str(SRS / script), *map(str, args)],
            check=True,
            cwd=ROOT,
            capture_output=True,
            text=True,
        )

    def test_cnvpytor_parser_retains_calls_and_soft_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            calls = tmp / "calls.txt"
            calls.write_text(
                "deletion\tchr1:100-1100\t1001\t0.50\t0.01\t0.02\t0.03\t0.04\t0.10\t0.00\t500000\n"
                "duplication\tchr2:200-2200\t2001\t1.60\t0.01\t0.02\t0.03\t0.04\t0.80\t0.70\t50000\n",
                encoding="utf-8",
            )
            genotypes = tmp / "gt.txt"
            genotypes.write_text(
                "chr1:100-1100\t0.50\t0.01\t0.02\t0.10\t0.00\t500000\t1\t2\t8\t0.30\t0.001\n",
                encoding="utf-8",
            )
            output = tmp / "out.tsv"
            self.run_script(
                "parse_cnvpytor_calls.py",
                "--calls", calls,
                "--genotypes", genotypes,
                "--sample", "S1",
                "--output", output,
            )
            rows = read_tsv(output)
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["SVTYPE"], "DEL")
            self.assertEqual(rows[0]["EVIDENCE_STATUS"], "RETAINED")
            self.assertEqual(rows[0]["CNVPYTOR_BAF_SHIFT"], "0.30")
            self.assertIn("HIGH_Q0_FRACTION", rows[1]["EVIDENCE_FLAGS"])
            self.assertIn("HIGH_REFERENCE_N_FRACTION", rows[1]["EVIDENCE_FLAGS"])
            self.assertIn("NEAR_LARGE_REFERENCE_GAP", rows[1]["EVIDENCE_FLAGS"])

    def test_gridss_reciprocal_breakends_are_one_event(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            vcf = tmp / "gridss.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                "chr1\t100\tbnd1\tN\tN]chr1:1000]\t100\tPASS\tMATEID=bnd2;EVENT=e1;AS=5;RAS=2\tSR:RP:VF\t4:8:10\n"
                "chr1\t1000\tbnd2\tN\t[chr1:100[N\t90\tPASS\tMATEID=bnd1;EVENT=e1;AS=5;RAS=2\tSR:RP:VF\t4:8:10\n",
                encoding="utf-8",
            )
            output = tmp / "events.tsv"
            self.run_script("gridss_vcf_to_events.py", "--vcf", vcf, "--output", output)
            rows = read_tsv(output)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["GRIDSS_EVENT_ID"], "e1")
            self.assertEqual(rows[0]["GRIDSS_RECORD_IDS"], "bnd1;bnd2")
            self.assertEqual(rows[0]["AS"], "5")
            self.assertEqual(rows[0]["SR"], "4")

    def test_mity_vaf_is_derived_from_allelic_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            vcf = tmp / "mity.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                "chrM\t100\t.\tA\tG\t100\tPASS\t.\tDP:RO:AO:VAF:tier:q\t1000:900:100:0.08:2:50\n",
                encoding="utf-8",
            )
            output = tmp / "mity.tsv"
            self.run_script(
                "parse_mity_vcf.py",
                "--vcf", vcf,
                "--sample", "S1",
                "--output", output,
            )
            rows = read_tsv(output)
            self.assertEqual(len(rows), 1)
            self.assertAlmostEqual(float(rows[0]["HETEROPLASMY_VAF"]), 0.08)
            self.assertEqual(rows[0]["MITY_TIER"], "2")
            self.assertEqual(rows[0]["MITY_Q"], "50")

    def test_srs_qc_uses_cumulative_mosdepth_distribution(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            flagstat = tmp / "flagstat.txt"
            flagstat.write_text(
                "1000 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "100 + 0 duplicates\n"
                "950 + 0 mapped (95.00% : N/A)\n"
                "900 + 0 properly paired (90.00% : N/A)\n",
                encoding="utf-8",
            )
            stats = tmp / "stats.txt"
            stats.write_text(
                "SN\taverage length:\t150\n"
                "SN\tinsert size average:\t350\n"
                "SN\tinsert size standard deviation:\t50\n"
                "SN\terror rate:\t0.001\n",
                encoding="utf-8",
            )
            summary = tmp / "mosdepth.summary.txt"
            summary.write_text(
                "chrom\tlength\tbases\tmean\tmin\tmax\n"
                "total\t1000\t30000\t30.0\t0\t80\n",
                encoding="utf-8",
            )
            dist = tmp / "mosdepth.global.dist.txt"
            dist.write_text(
                "total\t0\t1.0\n"
                "total\t10\t0.98\n"
                "total\t20\t0.90\n"
                "total\t29\t0.55\n"
                "total\t30\t0.49\n",
                encoding="utf-8",
            )
            output = tmp / "qc.tsv"
            self.run_script(
                "summarize_srs_qc.py",
                "--sample", "S1",
                "--flagstat", flagstat,
                "--stats", stats,
                "--mosdepth-summary", summary,
                "--mosdepth-global-dist", dist,
                "--output", output,
            )
            row = read_tsv(output)[0]
            self.assertEqual(row["MEDIAN_COVERAGE_FROM_GLOBAL_DIST"], "29")
            self.assertAlmostEqual(float(row["PERCENT_BASES_GE_10X"]), 98.0)
            self.assertAlmostEqual(float(row["PERCENT_BASES_GE_20X"]), 90.0)
            self.assertAlmostEqual(float(row["PERCENT_BASES_GE_30X"]), 49.0)

    def test_augment_keeps_master_and_emits_unmatched_cnv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            integrated = tmp / "integrated.tsv"
            integrated.write_text(
                "SV_ID\tCHROM\tSTART\tEND\tCHR2\tPOS2\tSVTYPE\tCALLERS\tGENES\n"
                "sv1\tchr1\t100\t1100\t.\t.\tDEL\tManta;Delly\tOPA1\n",
                encoding="utf-8",
            )
            cnv = tmp / "cnv.tsv"
            cnv.write_text(
                "SAMPLE\tCALLER\tSV_ID\tCHROM\tSTART\tEND\tSVTYPE\tSVLEN\tCNVPYTOR_LEVEL\tCNVPYTOR_EVAL1\tCNVPYTOR_Q0\tCNVPYTOR_PN\tCNVPYTOR_DG\tCNVPYTOR_BAF_SHIFT\tCNVPYTOR_BAF_PVALUE\tCNVPYTOR_HET_COUNT\tEVIDENCE_FLAGS\n"
                "S1\tCNVpytor\tc1\tchr1\t100\t1100\tDEL\t-1001\t0.5\t0.01\t0.1\t0.0\t500000\t0.3\t0.001\t8\t.\n"
                "S1\tCNVpytor\tc2\tchr2\t100\t5000\tDUP\t4901\t1.7\t0.01\t0.1\t0.0\t500000\t0.2\t0.01\t5\t.\n",
                encoding="utf-8",
            )
            ann = tmp / "ann.tsv"
            ann.write_text(
                "SV_ID\tGene_name\n"
                "c1\tOPA1\n"
                "c2\tMFN2\n",
                encoding="utf-8",
            )
            gridss = tmp / "gridss.tsv"
            gridss.write_text(
                "GRIDSS_EVENT_ID\tGRIDSS_RECORD_IDS\tCHROM1\tPOS1\tCHROM2\tPOS2\tORIENTATION\tQUAL_MAX\tFILTERS\tAS\tRAS\tCAS\tASSR\tASRP\tSR\tRP\tVF\n"
                "g1\tb1;b2\tchr1\t100\tchr1\t1100\t++\t100\tPASS\t5\t2\t.\t.\t.\t4\t8\t10\n",
                encoding="utf-8",
            )
            final = tmp / "final.tsv"
            cnv_only = tmp / "cnv_only.tsv"
            gridss_only = tmp / "gridss_only.tsv"
            self.run_script(
                "augment_srs_evidence.py",
                "--integrated", integrated,
                "--cnvpytor", cnv,
                "--cnvpytor-annotsv", ann,
                "--gridss-events", gridss,
                "--output", final,
                "--cnvpytor-only-output", cnv_only,
                "--gridss-only-output", gridss_only,
            )
            rows = read_tsv(final)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["SV_ID"], "sv1")
            self.assertEqual(rows[0]["CNVPYTOR_ID"], "c1")
            self.assertEqual(rows[0]["GRIDSS_EVENT_IDS"], "g1")
            self.assertIn("READ_DEPTH_CNV", rows[0]["SRS_EVIDENCE_MODALITIES"])
            self.assertIn("BAF_CONTEXT", rows[0]["SRS_EVIDENCE_MODALITIES"])
            unmatched = read_tsv(cnv_only)
            self.assertEqual(len(unmatched), 1)
            self.assertEqual(unmatched[0]["SV_ID"], "c2")
            self.assertEqual(unmatched[0]["ANNOTSV_GENES"], "MFN2")
            self.assertEqual(read_tsv(gridss_only), [])


if __name__ == "__main__":
    unittest.main()
