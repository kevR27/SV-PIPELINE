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
        self.assertEqual(config["survivor_caller_order"], "Manta,Delly")
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
