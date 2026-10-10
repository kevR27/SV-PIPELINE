import importlib.util
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]
SRS = ROOT / "snakemake_pipelines" / "srs_wgs_pipeline"
SCRIPTS = ROOT / "scripts"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


class SrsWgsPipelineTests(unittest.TestCase):
    def test_preflight_rejects_mismatched_contigs_and_empty_resources(self):
        script = module("preflight_srs", SCRIPTS / "preflight_srs_wgs.py")
        with tempfile.TemporaryDirectory() as tmp:
            bam_path = str(Path(tmp) / "sample.bam")
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": "chr1", "LN": 10000}]}
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                read = pysam.AlignedSegment()
                read.query_name = "paired_read"
                read.query_sequence = "A" * 50
                read.flag = 65
                read.reference_id = 0
                read.reference_start = 100
                read.mapping_quality = 60
                read.cigarstring = "50M"
                bam.write(read)
            pysam.index(bam_path)
            fai = Path(tmp) / "reference.fa.fai"
            for length, expected in [(10000, "PASS"), (9999, "FAIL")]:
                fai.write_text(f"chr1\t{length}\t6\t60\t61\n")
                checks = []
                script.check_bam_reference(bam_path, bam_path + ".bai", fai, checks)
                self.assertEqual(checks[-1]["status"], expected)
            empty_cache = Path(tmp) / "empty_cache"
            empty_cache.mkdir()
            checks = []
            script.record_file_check(str(empty_cache), "empty VEP cache", checks)
            self.assertEqual(checks[-1]["status"], "FAIL")

    def test_diagnostic_review_keeps_all_candidates_and_support_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "input.tsv", Path(tmp) / "output.tsv"
            with source.open("w") as handle:
                handle.write("SV_ID\tGENE\tTECHNICAL_SUPPORT\n")
                for n in range(104):
                    support = "READ_DEPTH_ONLY_CNV" if n == 0 else "SINGLE_CALLER"
                    if n == 103:
                        support = "MULTI_CALLER_PLUS_READ_DEPTH"
                    handle.write(f"sv{n}\tOTHER\t{support}\n")
            subprocess.run([sys.executable, str(SCRIPTS / "build_diagnostic_candidates.py"),
                            "--input", str(source), "--output", str(output)], check=True,
                           capture_output=True)
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 104)
            self.assertEqual(rows[0]["SV_ID"], "sv103")
            self.assertEqual(rows[-1]["SV_ID"], "sv0")

    def test_required_interfaces_are_wired(self):
        snake = (SRS / "Snakefile_SRS_WGS").read_text(encoding="utf-8")
        config = (SRS / "config_srs_wgs.yaml").read_text(encoding="utf-8")
        self.assertIn("reference/hon_hpo_seed_terms.tsv", config)
        self.assertIn("genes=rules.extract_sv_genes.output.genes", snake)
        self.assertIn("--hpo-seeds {input.hpo}", snake)
        self.assertIn("--edges {input.edges}", snake)
        self.assertIn("ANNOTSV_AUDIT_VCF", snake)
        self.assertIn(".annotsv.panel_only.tsv", snake)
        self.assertIn("_SV_VEP.panel_only.txt", snake)
        self.assertIn("{CANONICAL_FLAG}", snake)

    def test_missing_cnvpytor_qc_never_passes(self):
        script = module("combine_srs", SCRIPTS / "combine_sv_cnv_evidence.py")
        with tempfile.TemporaryDirectory() as tmp:
            calls = Path(tmp) / "calls.tsv"
            calls.write_text("deletion\tchr1:100-200\t101\t0.5\t.\t.\t.\t.\t.\t.\t1000\n", encoding="utf-8")
            rows = script.read_cnvpytor_calls(str(calls), 0.5, 0.5, 1e-4)
        self.assertEqual(rows[0]["status"], "REVIEW")
        self.assertIn("MISSING_Q0", rows[0]["flags"])
        self.assertIn("MISSING_PN", rows[0]["flags"])
        self.assertIn("MISSING_EVAL1", rows[0]["flags"])

    def test_gridss_rejects_failed_filter_and_low_qual(self):
        script = module("support_srs", SCRIPTS / "add_supporting_evidence.py")
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "gridss.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t100\tbad\tN\tN]chr1:200]\t99\tFAIL\tSVTYPE=BND\n"
                "chr1\t110\tlow\tN\tN]chr1:210]\t5\tPASS\tSVTYPE=BND\n"
                "chr1\t120\tgood\tN\tN]chr1:220]\t50\tPASS\tSVTYPE=BND\n",
                encoding="utf-8",
            )
            index = script.load_gridss(str(vcf), min_qual=10)
        self.assertEqual([item[3]["id"] for item in index["chr1"]], ["good"])

    def test_diagnostic_table_prioritizes_nuclear_mito(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "input.tsv", Path(tmp) / "output.tsv"
            source.write_text(
                "SV_ID\tGENE\tMITOCARTA_ENCODING\tPANEL_STATUS\tALLELE_RESEARCH_SCORE\n"
                "sv1\tOTHER\tNOT_MITOCARTA\tPANEL_GENE\t99\n"
                "sv2\tOPA1\tNUCLEAR_MITOCHONDRIAL_GENE\tPANEL_GENE\t1\n"
                "sv3\tMT-ND1\tMTDNA_ENCODED_GENE\tNON_PANEL\t100\n",
                encoding="utf-8",
            )
            subprocess.run([sys.executable, str(SCRIPTS / "build_diagnostic_candidates.py"), "--input", str(source), "--output", str(output)], check=True)
            lines = output.read_text(encoding="utf-8").splitlines()
        self.assertIn("NUCLEAR_MITOCHONDRIAL_PRIMARY", lines[1])
        self.assertIn("MTDNA_SECONDARY", lines[-1])


if __name__ == "__main__":
    unittest.main()
