"""Tests for separate Genomiser+CADD ranking and clinical-table integration."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from merge_genomiser_cadd_comparison import merge_patient, merge_cohort, check_same_inputs
from run_genomiser import build_java_command, validate_analysis


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class CaddComparisonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.base = self.root / "baseline.json"
        self.cadd = self.root / "cadd.json"
        common = {
            "sample_id": "14996",
            "small_vcf": "/data/14996/snv.vcf.gz",
            "sv_vcf": "/data/14996/sv.vcf.gz",
            "phenopacket_sha256": "matching-profile"
        }
        self.base.write_text(json.dumps({
            **common, "pathogenicity_sources": ["REMM", "SPLICE_AI"],
            "analysis_yaml_sha256": "baseline_hash"
        }))
        self.cadd.write_text(json.dumps({
            **common, "pathogenicity_sources": ["REMM", "SPLICE_AI", "CADD"],
            "analysis_yaml_sha256": "cadd_hash"
        }))

    def test_analysis_yaml_baseline_and_cadd_differ_only_in_sources(self):
        baseline = validate_analysis(ROOT / "config/genomiser/genome_analysis.yml")
        cadd = validate_analysis(ROOT / "config/genomiser/genome_analysis_cadd.yml")
        old_sources = baseline.pop("pathogenicitySources")
        new_sources = cadd.pop("pathogenicitySources")
        self.assertEqual(baseline, cadd)
        self.assertEqual(new_sources, old_sources + ["CADD"])

    def test_cadd_java_properties_are_explicit(self):
        args = SimpleNamespace(
            java="java", memory_gb=32, jar="/software/exomiser-cli.jar",
            data_dir="/data/exomiser", remm="/data/remm.tsv.gz",
            genome_version="2512", phenotype_version="2512",
            phenopacket="/data/phenopacket.json",
            analysis="/config/genome_analysis_cadd.yml",
            cadd_snv="/data/cadd/snvs.tsv.gz",
            cadd_indel="/data/cadd/indels.tsv.gz",
        )
        cmd = build_java_command(
            args, Path("/tmp/combined.vcf.gz"), Path("/tmp/options.yml")
        )
        self.assertIn("-Dexomiser.hg38.cadd-snv-path=/data/cadd/snvs.tsv.gz", cmd)
        self.assertIn("-Dexomiser.hg38.cadd-in-del-path=/data/cadd/indels.tsv.gz", cmd)

    def test_reject_different_patient_inputs(self):
        wrong = json.loads(self.cadd.read_text())
        wrong["phenopacket_sha256"] = "different-phenotype"
        self.cadd.write_text(json.dumps(wrong))
        with self.assertRaisesRegex(ValueError, "phenopacket_sha256"):
            check_same_inputs(self.base, self.cadd)

    def test_clinical_table_retains_baseline_and_adds_cadd(self):
        clinical = self.root / "clinical.tsv"
        genes = self.root / "cadd.genes.tsv"
        variants = self.root / "cadd.variants.tsv"
        output = self.root / "comparison.tsv"
        write_tsv(
            clinical,
            ["GENE", "CHROM", "START", "END", "SVTYPE",
             "GENOMISER_GENE_RANK", "FINAL_CLINICAL_CLASS"],
            [{"GENE": "OPA1", "CHROM": "chr3", "START": "100",
              "END": "200", "SVTYPE": "DEL", "GENOMISER_GENE_RANK": "8",
              "FINAL_CLINICAL_CLASS": "REVIEW"},
             {"GENE": "OTHER", "CHROM": "chr2", "START": "300",
              "END": "400", "SVTYPE": "DEL", "GENOMISER_GENE_RANK": "15",
              "FINAL_CLINICAL_CLASS": "REVIEW"}],
        )
        write_tsv(
            genes,
            ["#RANK", "GENE_SYMBOL", "MOI", "EXOMISER_GENE_COMBINED_SCORE"],
            [{"#RANK": "2", "GENE_SYMBOL": "OPA1", "MOI": "AD",
              "EXOMISER_GENE_COMBINED_SCORE": "0.8"}],
        )
        write_tsv(
            variants,
            ["GENE_SYMBOL", "CONTIG", "START", "END", "ALT", "VARIANT_SCORE"],
            [{"GENE_SYMBOL": "OPA1", "CONTIG": "3", "START": "100",
              "END": "200", "ALT": "<DEL>", "VARIANT_SCORE": "0.7"}],
        )
        merge_patient(clinical, genes, variants, self.base, self.cadd, output)
        with output.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["GENOMISER_GENE_RANK"], "8")
        self.assertEqual(rows[0]["GENOMISER_CADD_GENE_RANK"], "1")
        self.assertEqual(rows[0]["GENOMISER_CADD_RANK_IMPROVEMENT"], "7")
        self.assertEqual(rows[0]["FINAL_CLINICAL_CLASS"], "REVIEW")
        self.assertEqual(rows[0]["GENOMISER_CADD_SV_MATCH_STATUS"],
                         "EXACT_COORDINATE_CANDIDATE")
        self.assertEqual(rows[1]["GENOMISER_CADD_RANK_IMPROVEMENT"], ".")
        cohort = self.root / "cohort.tsv"
        merge_cohort(["14996=" + str(output)], cohort)
        with cohort.open(newline="") as handle:
            cohort_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(cohort_rows), 2)
        self.assertEqual(cohort_rows[0]["SAMPLE_ID"], "14996")


if __name__ == "__main__":
    unittest.main()
