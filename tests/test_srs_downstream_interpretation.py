"""Regression checks for SRS discovery retention and candidate interpretation."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
SRS = ROOT / "snakemake_pipelines" / "srs_wgs_pipeline"
SCRIPTS = ROOT / "scripts"


def load_script(filename):
    spec = importlib.util.spec_from_file_location(filename.removesuffix(".py"), SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SrsDownstreamInterpretationTests(unittest.TestCase):
    def test_single_caller_sv_union_remains_discovery_master(self):
        config = yaml.safe_load((SRS / "config_srs_wgs.yaml").read_text())
        snake = (SRS / "Snakefile_SRS_WGS").read_text()

        self.assertEqual(config["survivor_min_support"], 1)
        self.assertEqual(config["multi_caller_min_support"], 2)
        self.assertIn("SURVIVOR merge {output.list_file:q}", snake)
        self.assertIn("{SURVIVOR_MIN_SUPPORT}", snake)
        self.assertIn("combine_sv_and_cnv_calls.output.vcf", snake)
        self.assertIn('include: "srs_interpretation.smk"', snake)

    def test_srs_review_outputs_are_in_default_target(self):
        snake = (SRS / "Snakefile_SRS_WGS").read_text()
        module = (SRS / "srs_interpretation.smk").read_text()
        for suffix in (
            "_sv_gene_events.tsv.gz",
            "_large_sv_depth.tsv",
            "_large_sv_genes.tsv",
            "_snv_sv_candidates.tsv",
            "_gene_candidates.ranked.tsv",
            "_sv_gene_candidates.ranked.tsv",
            "_clinical_candidates.ranked.tsv",
        ):
            self.assertIn(suffix, snake)
        for rule in (
            "srs_interpret_sv_gene_events",
            "srs_large_sv_depth",
            "srs_snv_sv_candidates",
            "srs_finalize_candidate_ranking",
        ):
            self.assertIn(f"rule {rule}:", module)

    def test_snv_sv_pair_does_not_infer_trans_without_sv_phase(self):
        script = load_script("build_snv_sv_candidates.py")
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            vep = folder / "vep.tsv"
            phase = folder / "phased.vcf"
            sv = folder / "sv.tsv"
            output = folder / "pairs.tsv"

            vep.write_text(
                "#Uploaded_variation\tLocation\tAllele\tConsequence\tExisting_variation\tExtra\n"
                "chr1_100_A/G\t1:100\tG\tmissense_variant\t.\tSYMBOL=OPA1;IMPACT=MODERATE\n"
            )
            phase.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
                "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:PS\t0|1:42\n"
            )
            pd.DataFrame([{
                "GENE": "OPA1", "SV_ID": "sv01", "SVTYPE": "DEL",
                "CHROM": "chr1", "SV_GENE_EFFECT": "PARTIAL_GENE_DELETION",
                "POPULATION_STATUS": "UNKNOWN", "CALL_SUPPORT": "SINGLE_CALLER",
                "GENE_INHERITANCE_CLASS": "AUTOSOMAL_RECESSIVE",
                "GENE_MOI_SET": "AR",
                "INHERITANCE_MECHANISM_CLASS": "UNRESOLVED",
            }]).to_csv(sv, sep="\t", index=False)

            with patch.object(sys, "argv", [
                "build_snv_sv_candidates.py",
                "--vep", str(vep),
                "--phased-vcf", str(phase),
                "--sv-candidates", str(sv),
                "--output", str(output),
            ]):
                script.main()

            pairs = pd.read_csv(output, sep="\t", dtype=str)
            self.assertEqual(len(pairs), 1)
            self.assertEqual(pairs.iloc[0]["SMALL_PHASED"], "YES")
            self.assertEqual(pairs.iloc[0]["PHASE_RELATION"], "UNRESOLVED")
            self.assertEqual(
                pairs.iloc[0]["RECESSIVE_PAIR_STATUS"],
                "AR_SECOND_ALLELE_CANDIDATE_PHASE_UNRESOLVED",
            )

    def test_failed_or_uncalled_small_alleles_do_not_support_second_hit(self):
        script = load_script("build_snv_sv_candidates.py")
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            vep, phase, sv, output = [folder / name for name in
                                      ("vep.tsv", "phase.vcf", "sv.tsv", "pairs.tsv")]
            vep.write_text(
                "#Uploaded_variation\tLocation\tAllele\tConsequence\tExtra\n"
                + "".join(f"chr1_{pos}_A/G\t1:{pos}\tG\tmissense_variant\tSYMBOL=OPA1;IMPACT=MODERATE\n"
                          for pos in (100, 101, 102, 103))
                + "chr1_104_A/C,G\t1:104\tC\tmissense_variant\tSYMBOL=OPA1;IMPACT=MODERATE;ALLELE_NUM=1\n"
                + "chr1_104_A/C,G\t1:104\tG\tmissense_variant\tSYMBOL=OPA1;IMPACT=MODERATE;ALLELE_NUM=2\n"
            )
            phase.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr1\t100\t.\tA\tG\t10\tLowQual\t.\tGT:PS\t0|1:42\n"
                "chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT:PS\t0|0:42\n"
                "chr1\t102\t.\tA\tG\t50\tPASS\t.\tGT:PS\t.|.:42\n"
                "chr1\t103\t.\tA\tG\t50\tPASS\t.\tGT:PS\t0|1:42\n"
                "chr1\t104\t.\tA\tC,G\t50\tPASS\t.\tGT:PS\t0|1:42\n"
            )
            pd.DataFrame([{
                "GENE": "OPA1", "SV_ID": "sv01", "SVTYPE": "DEL", "CHROM": "chr1",
                "SV_GENE_EFFECT": "PARTIAL_GENE_DELETION", "POPULATION_STATUS": "UNKNOWN",
                "CALL_SUPPORT": "SINGLE_CALLER", "GENE_INHERITANCE_CLASS": "AUTOSOMAL_RECESSIVE",
                "GENE_MOI_SET": "AR", "LONGPHASE_GT": "1|0", "LONGPHASE_PS": "42",
                "LONGPHASE_PHASED": "YES",
            }]).to_csv(sv, sep="\t", index=False)
            with patch.object(sys, "argv", ["pairs", "--vep", str(vep), "--phased-vcf", str(phase),
                                           "--sv-candidates", str(sv), "--output", str(output)]):
                script.main()
            rows = pd.read_csv(output, sep="\t", dtype=str).set_index(["POS", "ALT"])
            for pos in ("100", "101", "102"):
                self.assertEqual(rows.loc[(pos, "G"), "RECESSIVE_PAIR_STATUS"], "SMALL_VARIANT_REQUIRES_REVIEW")
                self.assertEqual(rows.loc[(pos, "G"), "PHASE_RELATION"], "UNRESOLVED")
            self.assertEqual(rows.loc[("103", "G"), "RECESSIVE_PAIR_STATUS"], "AR_TRANS_SNV_SV_CANDIDATE")
            self.assertEqual(rows.loc[("104", "C"), "SMALL_CALL_STATUS"], "CALLED_ALT_PASS")
            self.assertEqual(rows.loc[("104", "G"), "SMALL_CALL_STATUS"], "ALT_NOT_IN_GENOTYPE")
            self.assertEqual(rows.loc[("104", "C"), "PHASE_RELATION"], "UNRESOLVED")

    def test_melt_support_requires_pass_and_called_alt(self):
        script = load_script("add_supporting_evidence.py")
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "melt.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr1\t100\tgood\tN\t<INS:ME:ALU>\t50\tPASS\tSVTYPE=INS\tGT\t0/1\n"
                "chr1\t101\tfailed\tN\t<INS:ME:ALU>\t10\tLowQual\tSVTYPE=INS\tGT\t0/1\n"
                "chr1\t102\treference\tN\t<INS:ME:ALU>\t50\tPASS\tSVTYPE=INS\tGT\t0/0\n"
                "chr1\t103\tmissing\tN\t<INS:ME:ALU>\t50\tPASS\tSVTYPE=INS\tGT\t./.\n"
            )
            records = script.load_melt(str(vcf))
            self.assertEqual([r["id"] for r in records["chr1"]], ["good"])

    def test_melt_single_receives_reference_list_not_zip(self):
        melt = load_script("run_melt_single.py")
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for filename in ("MELT.jar", "hg38.fa", "hg38.genes.bed", "ALU_MELT.zip", "sample.bam"):
                (folder / filename).write_bytes(b"test")
            output = folder / "family.vcf"
            observed = {}

            def fake_run(command, check):
                self.assertTrue(check)
                list_file = Path(command[command.index("-t") + 1])
                observed["reference_list"] = list_file.read_text()
                workdir = Path(command[command.index("-w") + 1])
                (workdir / "ALU.final_comp.vcf").write_text(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                )

            with patch.object(sys, "argv", [
                "run_melt_single.py",
                "--jar", str(folder / "MELT.jar"),
                "--bam", str(folder / "sample.bam"),
                "--reference", str(folder / "hg38.fa"),
                "--genes-bed", str(folder / "hg38.genes.bed"),
                "--transposon", str(folder / "ALU_MELT.zip"),
                "--workdir", str(folder / "work"),
                "--output", str(output),
            ]), patch.object(melt.subprocess, "run", side_effect=fake_run):
                self.assertEqual(melt.main(), 0)

            self.assertEqual(observed["reference_list"].strip(), str(folder / "ALU_MELT.zip"))
            self.assertTrue(output.is_file())


if __name__ == "__main__":
    unittest.main()
