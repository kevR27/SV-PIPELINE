"""Synthetic tests for patient phenotypes and Genomiser benchmarking."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_patient_profiles import load_profile, make_profile, terms
from benchmark_genomiser import chrom, gene_ranks, sv_match


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


class BenchmarkTests(unittest.TestCase):
    def test_chromosome_normalization(self):
        self.assertEqual(chrom("chrMT"), "MT")
        self.assertEqual(chrom("chrM"), "MT")
        self.assertEqual(chrom("chr3"), "3")

    def test_gene_ranks_deduplicates_inheritance_models(self):
        data = [{"GENE_SYMBOL": "OPA1"}, {"GENE_SYMBOL": "OPA1"}, {"GENE_SYMBOL": "ACO2"}]
        self.assertEqual(gene_ranks(data, "GENE_SYMBOL"), {"OPA1": 1, "ACO2": 2})

    def test_sv_event_coordinate_match(self):
        truth = {"gene": "OPA1", "chrom": "chr3", "pos": "1000", "end": "2000", "svtype": "DEL"}
        custom = {"GENE": "OPA1", "CHROM": "chr3", "START": "1005",
                  "END": "1990", "SVTYPE": "DEL"}
        genomiser = {"GENE_SYMBOL": "OPA1", "CONTIG": "3", "START": "1000",
                     "END": "2000", "ALT": "<DEL>"}
        self.assertTrue(sv_match(custom, truth, "custom", 20))
        self.assertTrue(sv_match(genomiser, truth, "genomiser", 20))
        self.assertFalse(sv_match(genomiser, {**truth, "svtype": "DUP"}, "genomiser", 20))


if __name__ == "__main__":
    unittest.main()
