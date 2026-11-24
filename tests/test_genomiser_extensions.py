"""Synthetic tests for patient phenotypes and Genomiser benchmarking."""
import csv
import json
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_patient_profiles import load_profile, make_profile, terms, main as prepare_profiles_main
from genomiser_requirements import require_genomiser
from run_genomiser import build_java_command, validate_analysis
from benchmark_genomiser import (chromosome, sequential_gene_ranks,
                                 custom_global_gene_ranks, genomiser_global_gene_ranks,
                                 sv_matches)


class PatientProfileTests(unittest.TestCase):
    def base(self):
        return {
            "sample_id": "TEST01", "hpo_present": "HP:0000648;HP:0007663",
            "hpo_absent": "HP:0000543", "onset_min_years": "12",
            "onset_max_years": "16", "onset_status": "estimated",
            "diagnostic_hypothesis": "suspected dominant optic atrophy",
            "hypothesis_status": "suspected", "sex": "UNKNOWN_SEX"
        }

    def test_valid_profile_keeps_hypothesis_separate(self):
        packet, clinical = make_profile(self.base(), "TEST01", {})
        self.assertEqual([x["type"]["id"] for x in packet["phenotypicFeatures"]],
                         ["HP:0000648", "HP:0007663"])
        self.assertEqual(clinical["hpo_absent"], ["HP:0000543"])
        self.assertEqual(clinical["onset_min_years"], 12)
        self.assertNotIn("diseases", packet)

    def test_shared_hon_seeds_combine_with_observed_hpo(self):
        case = self.base()
        case["hpo_present"] += ";HP:0001250"
        seeds = [
            "HP:0000648", "HP:0007663", "HP:0000543",
            "HP:0000529", "HP:0000603", "HP:0000576"
        ]
        packet, context = make_profile(case, "TEST01", {}, hon_seed_terms=seeds)
        scored = {feature["type"]["id"] for feature in packet["phenotypicFeatures"]}
        self.assertEqual(
            scored,
            {"HP:0000648", "HP:0007663", "HP:0001250",
             "HP:0000529", "HP:0000603", "HP:0000576"}
        )
        self.assertEqual(context["hpo_present"],
                         ["HP:0000648", "HP:0001250", "HP:0007663"])
        self.assertEqual(context["hpo_hon_reference_excluded_as_absent"],
                         ["HP:0000543"])
        self.assertEqual(context["hpo_scoring_term_sources"]["HP:0000603"],
                         "SHARED_HON_REFERENCE")
        self.assertEqual(context["hpo_scoring_term_sources"]["HP:0001250"],
                         "PATIENT_OBSERVED")
        self.assertEqual(context["hpo_scoring_terms"], sorted(scored))

    def test_no_hon_profile_imputed_to_non_oa_control(self):
        case = self.base()
        case["hpo_present"] = "HP:0001250"
        case["hpo_absent"] = ""
        packet, context = make_profile(
            case, "TEST01", {}, hon_seed_terms=["HP:0000648", "HP:0007663"]
        )
        self.assertEqual([f["type"]["id"] for f in packet["phenotypicFeatures"]],
                         ["HP:0001250"])
        self.assertFalse(context["hpo_hon_seed_profile_applied"])

    def test_observed_absence_overrides_shared_seed(self):
        case = self.base()
        case["hpo_absent"] = "HP:0000603"
        packet, context = make_profile(
            case, "TEST01", {}, hon_seed_terms=["HP:0000648", "HP:0000603"]
        )
        self.assertNotIn("HP:0000603", context["hpo_scoring_terms"])
        self.assertEqual(context["hpo_hon_reference_excluded_as_absent"],
                         ["HP:0000603"])
        self.assertNotIn("HP:0000603",
                         [f["type"]["id"] for f in packet["phenotypicFeatures"]])

    def test_profile_files_use_same_union_for_genomiser_and_custom_ranking(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            clinical = root / "patients.tsv"
            seeds = root / "seeds.tsv"
            packet_file = root / "patient.phenopacket.json"
            hpo_file = root / "patient.patient_hpo.tsv"
            context_file = root / "patient.clinical_context.json"
            clinical.write_text(
                "sample_id\thpo_present\thpo_absent\n"
                "TEST01\tHP:0000648;HP:0001250\tHP:0000603\n",
                encoding="utf-8"
            )
            seeds.write_text(
                "hpo_id\thpo_label\thon_seed_role\n"
                "HP:0000648\tOptic atrophy\tCORE_OCULAR_HON\n"
                "HP:0007663\tReduced visual acuity\tCORE_OCULAR_HON\n"
                "HP:0000603\tCentral scotoma\tCORE_OCULAR_HON\n",
                encoding="utf-8"
            )
            arguments = [
                "prepare_patient_profiles.py",
                "--profiles", str(clinical),
                "--sample", "TEST01",
                "--hpo-seeds", str(seeds),
                "--include-shared-hon-seeds",
                "--phenopacket", str(packet_file),
                "--hpo-output", str(hpo_file),
                "--clinical-output", str(context_file),
            ]
            with patch.object(sys, "argv", arguments):
                prepare_profiles_main()
            packet = json.loads(packet_file.read_text())
            context = json.loads(context_file.read_text())
            with hpo_file.open(newline="") as handle:
                hpo_rows = list(csv.DictReader(handle, delimiter="\t"))
            packet_terms = {row["type"]["id"] for row in packet["phenotypicFeatures"]}
            hpo_terms = {row["hpo_id"] for row in hpo_rows}
            self.assertEqual(packet_terms, hpo_terms)
            self.assertEqual(packet_terms, {"HP:0000648", "HP:0001250", "HP:0007663"})
            self.assertEqual(context["hpo_present"], ["HP:0000648", "HP:0001250"])
            self.assertEqual(context["hpo_hon_reference_added"], ["HP:0007663"])
            self.assertEqual(context["hpo_hon_reference_excluded_as_absent"],
                             ["HP:0000603"])

    def test_default_does_not_impute_shared_hon_terms(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            clinical = root / "patients.tsv"
            seeds = root / "seeds.tsv"
            packet_file = root / "patient.phenopacket.json"
            hpo_file = root / "patient.patient_hpo.tsv"
            context_file = root / "patient.clinical_context.json"
            clinical.write_text(
                "sample_id\\thpo_present\\thpo_absent\\n"
                "TEST01\\tHP:0000648\\t\\n",
                encoding="utf-8",
            )
            seeds.write_text(
                "hpo_id\\thpo_label\\thon_seed_role\\n"
                "HP:0000648\\tOptic atrophy\\tCORE_OCULAR_HON\\n"
                "HP:0000603\\tCentral scotoma\\tCORE_OCULAR_HON\\n",
                encoding="utf-8",
            )
            args = [
                "prepare_patient_profiles.py", "--profiles", str(clinical),
                "--sample", "TEST01", "--hpo-seeds", str(seeds),
                "--phenopacket", str(packet_file),
                "--hpo-output", str(hpo_file),
                "--clinical-output", str(context_file),
            ]
            with patch.object(sys, "argv", args):
                prepare_profiles_main()
            packet = json.loads(packet_file.read_text())
            context = json.loads(context_file.read_text())
            self.assertEqual(
                [f["type"]["id"] for f in packet["phenotypicFeatures"]],
                ["HP:0000648"],
            )
            self.assertEqual(context["hpo_scoring_terms"], ["HP:0000648"])
            self.assertEqual(context["hpo_hon_reference_added"], [])
            self.assertEqual(context["hpo_absent"], [])

    def test_invalid_hpo_is_rejected(self):
        with self.assertRaises(ValueError):
            terms("HP:0000648;HELLO")

    def test_contradictory_hpo_is_rejected(self):
        case = self.base()
        case["hpo_absent"] = "HP:0000648"
        with self.assertRaises(ValueError):
            make_profile(case, "TEST01", {})

    def test_invalid_onset_is_rejected(self):
        case = self.base()
        case["onset_min_years"] = "17"
        with self.assertRaises(ValueError):
            make_profile(case, "TEST01", {})

    def test_duplicate_patient_ids_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "patients.tsv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["sample_id", "hpo_present"], delimiter="\t")
                writer.writeheader()
                writer.writerows([{"sample_id": "TEST01", "hpo_present": "HP:0000648"}] * 2)
            with self.assertRaises(ValueError):
                load_profile(path, "TEST01")


class VcfOnlyWorkflowTests(unittest.TestCase):
    def test_bam_free_entrypoint_contains_only_genomiser_rules(self):
        snakefile = (ROOT / "snakemake_pipelines" / "genomiser_vcf_only" /
                     "Snakefile").read_text(encoding="utf-8")
        self.assertIn('include: "../common/patient_profiles.smk"', snakefile)
        self.assertIn('include: "../common/genomiser.smk"', snakefile)
        self.assertNotIn("bam=lambda", snakefile)
        self.assertNotIn("bai=lambda", snakefile)
        self.assertNotIn("rule clair3_snp_calling:", snakefile)
        self.assertNotIn("rule jasmine_merge:", snakefile)
        self.assertNotIn("Snakefile_LRS_update", snakefile)
        self.assertIn('VCF_INPUTS[wc.sample]["small_vcf"]', snakefile)
        self.assertIn('VCF_INPUTS[wc.sample]["sv_vcf"]', snakefile)

    def test_vcf_only_config_has_no_bams(self):
        import yaml
        path = (ROOT / "snakemake_pipelines" / "genomiser_vcf_only" /
                "config.yaml")
        cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
        self.assertNotIn("samples", cfg)
        self.assertNotIn("bam", cfg)
        self.assertTrue(cfg["vcf_samples"])
        for sample, inputs in cfg["vcf_samples"].items():
            self.assertTrue(str(sample))
            self.assertTrue(inputs["small_vcf"].endswith(".vcf.gz"))
            self.assertTrue(inputs["sv_vcf"].endswith(".vcf.gz"))
            self.assertNotEqual(inputs["small_vcf"], inputs["sv_vcf"])
        self.assertEqual(cfg["genomiser_enabled"], True)
        self.assertTrue(cfg["patient_profiles_tsv"])
        self.assertTrue(cfg["genomiser_analysis_yaml"])


class BenchmarkTests(unittest.TestCase):
    def test_chromosome_normalization(self):
        self.assertEqual(chromosome("chrMT"), "MT")
        self.assertEqual(chromosome("chrM"), "MT")
        self.assertEqual(chromosome("chr3"), "3")

    def test_gene_ranks_deduplicates_inheritance_models(self):
        data = [{"GENE_SYMBOL": "OPA1"}, {"GENE_SYMBOL": "OPA1"}, {"GENE_SYMBOL": "ACO2"}]
        self.assertEqual(sequential_gene_ranks(data, "GENE_SYMBOL"), {"OPA1": 1, "ACO2": 2})

    def test_custom_gene_rank_is_global_not_panel_sorted(self):
        genes = [
            {"GENE": "LOW", "PANEL_STATUS": "PANEL", "FINAL_GENE_RELEVANCE_TIER": "LIMITED",
             "BEST_EVENT_FINAL_RANK_GLOBAL": "20", "FINAL_GENE_RELEVANCE_DISPLAY_SCORE": "0"},
            {"GENE": "HIGH", "PANEL_STATUS": "NONPANEL", "FINAL_GENE_RELEVANCE_TIER": "HIGH",
             "BEST_EVENT_FINAL_RANK_GLOBAL": "1", "FINAL_GENE_RELEVANCE_DISPLAY_SCORE": "9"}
        ]
        self.assertEqual(custom_global_gene_ranks(genes, set(genes[0])),
                         {"HIGH": 1, "LOW": 2})

    def test_genomiser_genes_ranked_by_model_order(self):
        rows = [{"GENE_SYMBOL": "OPA1", "#RANK": "2"},
                {"GENE_SYMBOL": "ACO2", "#RANK": "1"}]
        self.assertEqual(genomiser_global_gene_ranks(rows, set(rows[0])),
                         {"ACO2": 1, "OPA1": 2})

    def test_sv_event_coordinate_match(self):
        truth = {"gene": "OPA1", "chrom": "chr3", "pos": "1000", "end": "2000", "svtype": "DEL"}
        custom = {"GENE": "OPA1", "CHROM": "chr3", "START": "1005",
                  "END": "1990", "SVTYPE": "DEL"}
        genomiser = {"GENE_SYMBOL": "OPA1", "CONTIG": "3", "START": "1000",
                     "END": "2000", "ALT": "<DEL>"}
        self.assertTrue(sv_matches(custom, truth, "custom", 20))
        self.assertTrue(sv_matches(genomiser, truth, "genomiser", 20))
        self.assertFalse(sv_matches(genomiser, {**truth, "svtype": "DUP"}, "genomiser", 20))


class EditableGenomeAnalysisTests(unittest.TestCase):
    def baseline(self):
        return ROOT / "config" / "genomiser" / "genome_analysis.yml"

    def test_baseline_is_valid_genome_analysis(self):
        config = validate_analysis(self.baseline())
        self.assertEqual(config["analysisMode"], "PASS_ONLY")
        self.assertIn("REMM", config["pathogenicitySources"])
        self.assertIn("SPLICE_AI", config["pathogenicitySources"])

    def test_reject_exome_like_interval_filter(self):
        import yaml
        analysis = validate_analysis(self.baseline())
        analysis["steps"].append({"genePanelFilter": {"geneSymbols": ["OPA1"]}})
        with tempfile.TemporaryDirectory() as directory:
            filepath = Path(directory) / "analysis.yml"
            filepath.write_text(yaml.safe_dump(analysis), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "genePanelFilter"):
                validate_analysis(filepath)

    def test_reject_intron_removal(self):
        import yaml
        analysis = validate_analysis(self.baseline())
        analysis["steps"].append({
            "variantEffectFilter": {"remove": ["CODING_TRANSCRIPT_INTRON_VARIANT"]}
        })
        with tempfile.TemporaryDirectory() as directory:
            filepath = Path(directory) / "analysis.yml"
            filepath.write_text(yaml.safe_dump(analysis), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "intronic/regulatory"):
                validate_analysis(filepath)

    def test_uses_custom_analysis_not_preset(self):
        args = SimpleNamespace(
            java="java",
            memory_gb=32,
            jar="/soft/exomiser-cli-15.1.0.jar",
            data_dir="/data/exomiser",
            genome_version="2512",
            phenotype_version="2512",
            remm="/data/ReMM.v0.4.hg38.tsv.gz",
            phenopacket="/private/patient.phenopacket.json",
            analysis=str(self.baseline()),
        )
        command = build_java_command(
            args, Path("/data/one.genomiser.input.vcf.gz"),
            Path("/data/output-options.yml")
        )
        self.assertEqual(command[command.index("--analysis") + 1],
                         str(self.baseline().resolve()))
        self.assertNotIn("--preset", command)
        self.assertEqual(command[command.index("--assembly") + 1], "hg38")


class MandatoryGenomiserTests(unittest.TestCase):
    def settings(self):
        return {
            "genomiser_enabled": True,
            "patient_profiles_tsv": "/private/patients.tsv",
            "genomiser_cli_jar": "/software/exomiser-cli.jar",
            "genomiser_data_dir": "/data/exomiser",
            "genomiser_remm_hg38": "/data/remm/ReMM.v0.4.hg38.tsv.gz",
            "genomiser_memory_gb": 32,
            "genomiser_analysis_yaml": "/software/config/genomiser/genome_analysis.yml",
        }

    def test_mandatory_enabled_config(self):
        self.assertTrue(require_genomiser(self.settings(), ["21930"]))

    def test_explicit_disable_is_rejected(self):
        settings = self.settings()
        settings["genomiser_enabled"] = False
        with self.assertRaisesRegex(ValueError, "required"):
            require_genomiser(settings, ["21930"])

    def test_cli_string_true_is_accepted(self):
        settings = self.settings()
        settings["genomiser_enabled"] = "true"
        self.assertTrue(require_genomiser(settings, ["21930"]))

    def test_cli_string_false_is_rejected(self):
        settings = self.settings()
        settings["genomiser_enabled"] = "false"
        with self.assertRaisesRegex(ValueError, "required"):
            require_genomiser(settings, ["21930"])

    def test_missing_patient_profiles_is_rejected(self):
        settings = self.settings()
        settings["patient_profiles_tsv"] = None
        with self.assertRaisesRegex(ValueError, "patient_profiles_tsv"):
            require_genomiser(settings, ["21930"])

    def test_missing_remm_is_rejected(self):
        settings = self.settings()
        settings["genomiser_remm_hg38"] = ""
        with self.assertRaisesRegex(ValueError, "genomiser_remm_hg38"):
            require_genomiser(settings, ["21930"])

    def test_no_samples_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "at least one sample"):
            require_genomiser(self.settings(), [])

    def test_patient_hpo_override_is_rejected(self):
        settings = self.settings()
        settings["patient_hpo_files"] = {"21930": "/private/other.hpo"}
        with self.assertRaisesRegex(ValueError, "do not also"):
            require_genomiser(settings, ["21930"])


if __name__ == "__main__":
    unittest.main()
