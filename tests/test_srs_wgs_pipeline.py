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
            # A separately supplied BAI must also be discoverable by callers.
            fai.write_text("chr1\t10000\t6\t60\t61\n")
            custom_index = Path(tmp) / "elsewhere.bai"
            Path(bam_path + ".bai").rename(custom_index)
            checks = []
            script.check_bam_reference(bam_path, str(custom_index), fai, checks)
            self.assertEqual(checks[-1]["status"], "FAIL")
            Path(bam_path + ".bai").symlink_to(custom_index)
            checks = []
            script.check_bam_reference(bam_path, str(custom_index), fai, checks)
            self.assertEqual(checks[-1]["status"], "PASS")
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
        self.assertIn("--hpo-seeds {input.hpo:q}", snake)
        self.assertIn("--edges {input.edges:q}", snake)
        self.assertIn("ANNOTSV_AUDIT_VCF", snake)
        self.assertIn(".annotsv.panel_only.tsv", snake)
        self.assertIn("_SV_VEP.panel_only.txt", snake)
        self.assertIn("{CANONICAL_FLAG}", snake)

    def test_srs_sv_vep_uses_supported_buffer_size(self):
        snake = (SRS / "Snakefile_SRS_WGS").read_text(encoding="utf-8")
        rule = snake.split("rule vep_sv_annotation:", 1)[1].split("\nrule ", 1)[0]
        self.assertIn("--buffer_size 100", rule)
        self.assertNotIn("--batch_size", rule)
        self.assertIn("--max_sv_size 250000000", rule)
        self.assertIn("--fork {threads}", rule)

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

    def test_cnvpytor_checks_hg38_without_other_genomes(self):
        script = (SCRIPTS / "run_cnvpytor_hg38.py").read_text()
        snake = (SRS / "Snakefile_SRS_WGS").read_text()
        rule = snake.split("rule cnvpytor_cnv_calling:", 1)[1].split("\nrule ", 1)[0]

        self.assertIn('Genome.reference_genomes = {"hg38": hg38}', script)
        self.assertIn('for key in ("gc_file", "mask_file")', script)
        self.assertIn("cnvpytor_main()", script)
        self.assertIn('script=SCRIPTS + "/run_cnvpytor_hg38.py"', rule)
        self.assertNotIn("\n        cnvpytor -root ", rule)
        self.assertIn('test -s {output.pytor:q}', rule)
        self.assertIn('cat {output.info:q}', rule)
        self.assertIn("Usingreferencegenome:hg38[GC:yes,mask:yes]", rule)
        self.assertIn('logs/cnvpytor/{sample}.cnvpytor.log', rule)

    def test_deepvariant_optional_gpu_configuration(self):
        config = (SRS / "config_srs_wgs.yaml").read_text(encoding="utf-8")
        snake = (SRS / "Snakefile_SRS_WGS").read_text(encoding="utf-8")
        rule = snake.split("rule deepvariant_snv_indel_calling:", 1)[1].split("\nrule ", 1)[0]

        self.assertIn("deepvariant_use_gpu: false", config)
        self.assertIn("deepvariant_gpu_image: \"docker://google/deepvariant:1.10.0-gpu\"", config)
        self.assertIn('gpu_flag="--nv" if DEEPVARIANT_USE_GPU else ""', rule)
        self.assertIn("gpus=1 if DEEPVARIANT_USE_GPU else 0", rule)
        self.assertIn("image=DEEPVARIANT_ACTIVE_IMAGE", rule)
        self.assertIn("apptainer exec {params.gpu_flag}", rule)
        self.assertIn("--model_type=WGS", rule)

    def test_manta_run_script_uses_python2_on_noexec_filesystem(self):
        snake = (SRS / "Snakefile_SRS_WGS").read_text()
        manta_rule = snake.split("rule manta_sv_calling:", 1)[1].split("\nrule ", 1)[0]
        self.assertIn("python2 {params.rundir:q}/runWorkflow.py -m local -j {threads}", manta_rule)
        self.assertIn('python2 "$(command -v configManta.py)"', manta_rule)
        self.assertIn("command -v python2", manta_rule)
        self.assertIn("for tool in samtools bgzip tabix", manta_rule)
        self.assertIn("convertInversion.py", manta_rule)
        manta_env = (ROOT / "envs" / "manta.yaml").read_text()
        self.assertIn("  - python=2.7", manta_env)
        self.assertIn("  - bioconda::samtools=1.9", manta_env)
        self.assertNotIn("htslib=1.20", manta_env)
        self.assertNotIn("samtools=1.20", manta_env)
        self.assertNotIn("\n        {params.rundir:q}/runWorkflow.py", manta_rule)

    def test_gridss_snappy_native_library_uses_home_tmp(self):
        snake = (ROOT / "snakemake_pipelines" / "srs_wgs_pipeline" / "Snakefile_SRS_WGS").read_text()
        rule = snake.split("rule gridss_sv_calling:", 1)[1].split("\nrule ", 1)[0]
        self.assertIn("org.xerial.snappy.tempdir=$HOME/.cache/gridss-native", rule)
        self.assertIn("java.io.tmpdir=$HOME/.cache/gridss-native", rule)
        self.assertIn("mkdir -p \"$HOME/.cache/gridss-native\"", rule)
        self.assertIn('include: "srs_interpretation.smk"', snake)

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
