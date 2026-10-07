import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LRS = ROOT / "snakemake_pipelines" / "lrs"
SRS = ROOT / "snakemake_pipelines" / "srs_wgs_pipeline"


class PipelineParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lrs = (LRS / "Snakefile_LRS_update").read_text(encoding="utf-8")
        cls.lrs_post = (LRS / "Snakefile_LRS_postprocess").read_text(
            encoding="utf-8"
        )
        cls.srs = (SRS / "Snakefile_SRS_WGS").read_text(encoding="utf-8")

    def test_shared_interpretation_scripts_are_used(self):
        shared_scripts = [
            "filter_annotsv_panel.py",
            "extract_annotsv_genes.py",
            "extract_nonpanel_genes.py",
            "monarch_human_gene_phenotypes.py",
            "rank_sv_gene_candidates.py",
            "build_integrated_sv_gene_tsv.py",
        ]
        for script in shared_scripts:
            with self.subTest(script=script):
                self.assertIn(script, self.lrs)
                self.assertIn(script, self.srs)

    def test_shared_assessment_and_audit_layers_are_included(self):
        for include in ["allele_assessment.smk", "annotsv_evidence_audit.smk"]:
            with self.subTest(include=include):
                self.assertIn(include, self.lrs)
                self.assertIn(include, self.srs)

    def test_platform_specific_layers_remain_separate(self):
        self.assertIn("modkit_pileup", self.lrs)
        self.assertNotIn("modkit_pileup", self.srs)
        self.assertIn("gridss_sv_calling", self.srs)
        self.assertNotIn("gridss_sv_calling", self.lrs)
        self.assertIn("cnvpytor_cnv_calling", self.srs)
        self.assertNotIn("cnvpytor_cnv_calling", self.lrs)
        self.assertIn("needlr_annotation", self.lrs)
        self.assertNotIn("needlr_annotation", self.srs)

    def test_nuclear_mito_is_primary_and_mtdna_is_secondary(self):
        lrs_config = (LRS / "config_lrs.yaml").read_text(encoding="utf-8")
        srs_config = (SRS / "config_srs_wgs.yaml").read_text(encoding="utf-8")
        mito_ranker = (ROOT / "scripts" / "rank_mitochondrial_genes.py").read_text(
            encoding="utf-8"
        )

        self.assertIn('analysis_primary_focus: "nuclear_mitochondrial_genes"', lrs_config)
        self.assertIn('mtdna_analysis_role: "secondary"', lrs_config)
        self.assertIn("mitocarta_enabled: true", srs_config)
        self.assertIn("PRIMARY_NUCLEAR_MITOCHONDRIAL", mito_ranker)
        self.assertIn("SECONDARY_MTDNA", mito_ranker)
        self.assertIn("rank_mitochondrial_genes.py", self.lrs_post)
        self.assertIn("rank_mitochondrial_genes.py", self.srs)

    def test_lrs_mitochondrial_output_places_nuclear_genes_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            input_table = tmp / "input.tsv"
            gene_bed = tmp / "genes.bed"
            panel = tmp / "panel.txt"
            output = tmp / "output.tsv"

            input_table.write_text(
                "SV_ID\tGENE\tMITOCARTA_ENCODING\tSV_GENE_RELATIONSHIP\t"
                "GENE_RELEVANCE_TIER\tEVENT_GENE_RELEVANCE_SCORE\n"
                "sv_mt\tMT-ND1\tMTDNA_ENCODED_GENE\t"
                "WHOLE_GENE_DOSAGE_CONTEXT\tHIGH\t99\n"
                "sv_nuclear\tOPA1\tNUCLEAR_MITOCHONDRIAL_GENE\t"
                "PARTIAL_GENE_OVERLAP\tMODERATE\t1\n",
                encoding="utf-8",
            )
            gene_bed.write_text("chrM\t3306\t4262\tMT-ND1\n", encoding="utf-8")
            panel.write_text("OPA1\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "rank_mitochondrial_genes.py"),
                    "--input", str(input_table),
                    "--gene-bed", str(gene_bed),
                    "--panel", str(panel),
                    "--output", str(output),
                ],
                check=True,
                env={"PYTHONPATH": str(ROOT / "scripts")},
            )

            rows = output.read_text(encoding="utf-8").splitlines()

        self.assertIn("PRIMARY_NUCLEAR_MITOCHONDRIAL", rows[1])
        self.assertIn("SECONDARY_MTDNA", rows[2])

    def test_lrs_chrM_summary_states_its_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            vcf = tmp / "mt.vcf.gz"
            output = tmp / "mt.tsv"

            with gzip.open(vcf, "wt", encoding="utf-8") as handle:
                handle.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                    "chrM\t3460\tmt1\tA\tG\t50\tPASS\tDP=30\t"
                    "GT:DP:AF\t0/1:30:0.25\n"
                )

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "summarize_secondary_mtdna_vcf.py"),
                    "--vcf", str(vcf),
                    "--source", "Clair3 long-read candidate review",
                    "--output", str(output),
                ],
                check=True,
            )

            result = output.read_text(encoding="utf-8")

        self.assertIn("SECONDARY_MTDNA_REVIEW", result)
        self.assertIn("not a validated heteroplasmy", result)


if __name__ == "__main__":
    unittest.main()
