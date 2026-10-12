"""Synthetic regressions for the separate Exomiser/HON comparison."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from compare_exomiser_hon import backend_resource_paths, check_backend_files, exomiser_by_gene, extract_exomiser, hpo_ids
from sv_evidence_common import hon_hpo_seed_path


def table(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", restval=".")
        writer.writeheader()
        writer.writerows(rows)


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ExomiserComparisonTests(unittest.TestCase):
    def test_one_explicit_seed_path_and_older_configs(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "different_hon_profile.tsv"
            expected = str(path.resolve())
            self.assertEqual(hon_hpo_seed_path({"hon_hpo_seed_terms": str(path)}, ROOT), expected)
            self.assertEqual(hon_hpo_seed_path({"discovery_hpo_file": str(path)}, ROOT), expected)
            self.assertEqual(hon_hpo_seed_path({"hon_hpo_seed_terms": str(path), "discovery_hpo_file": str(path)}, ROOT), expected)
        self.assertEqual(hon_hpo_seed_path({}, ROOT), str(ROOT / "reference/hon_hpo_seed_terms.tsv"))

    def test_conflicting_or_invalid_seed_paths_stop_before_scoring(self):
        with self.assertRaisesRegex(ValueError, "different files"):
            hon_hpo_seed_path({"hon_hpo_seed_terms": "a.tsv", "discovery_hpo_file": "b.tsv"}, ROOT)
        for raw in ["", " ", True, ["a.tsv"]]:
            with self.subTest(raw=raw), self.assertRaisesRegex(ValueError, "nonempty path"):
                hon_hpo_seed_path({"hon_hpo_seed_terms": raw}, ROOT)

    def test_custom_seed_file_changes_the_saved_query(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            seed = path / "different_hon_profile.tsv"
            table(seed, [{"hpo_id": "HP:0000648", "hpo_label": "Optic atrophy", "hon_seed_role": "CORE_OCULAR_HON"},
                         {"hpo_id": "HP:0001251", "hpo_label": "Ataxia", "hon_seed_role": "MITO_SYNDROMIC_CONTEXT"}])
            resolved = hon_hpo_seed_path({"hon_hpo_seed_terms": str(seed)}, ROOT)
            self.command("prepare-hpo", "--enabled", "--hpo-seeds", resolved,
                         "--output", path / "query.tsv", "--terms-output", path / "query.txt")
            self.assertEqual((path / "query.txt").read_text(), "HP:0000648")
            self.assertEqual(read(path / "query.tsv")[0]["query_source"], str(seed))

    def command(self, *arguments, success=True):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/compare_exomiser_hon.py"), *map(str, arguments)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result

    def prepare(self, path, *extra):
        self.command("prepare-hpo", "--hpo-seeds", ROOT / "reference/hon_hpo_seed_terms.tsv",
                     "--output", path / "query.tsv", "--terms-output", path / "query.txt", *extra)

    def fixtures(self, path):
        self.prepare(path, "--enabled", "--shared-hpo-json", '"HP:0000648"')
        table(path / "rank.tsv", [
            {"gene": "A", "PANEL_STATUS": "PANEL_GENE", "GENE_RANK_GLOBAL": "1", "GENE_RELEVANCE_TIER": "HIGH", "GENE_RELEVANCE_DISPLAY_SCORE": "9.75", "HON_SEMANTIC_SCORE_0_10": "9.5"},
            {"gene": "B", "PANEL_STATUS": "NONPANEL_GENE", "GENE_RANK_GLOBAL": "2", "GENE_RELEVANCE_TIER": "MODERATE", "GENE_RELEVANCE_DISPLAY_SCORE": "6.25", "HON_SEMANTIC_SCORE_0_10": "5.0"},
            {"gene": "C", "PANEL_STATUS": "NONPANEL_GENE", "GENE_RANK_GLOBAL": "3", "GENE_RELEVANCE_TIER": "LIMITED", "GENE_RELEVANCE_DISPLAY_SCORE": "0.0", "HON_SEMANTIC_SCORE_0_10": "0.0"},
        ])
        table(path / "pheno.tsv", [{"gene_symbol": "A", "hpo_id": "HP:0000648"},
                                    {"gene_symbol": "B", "hpo_id": "HP:0007663"},
                                    {"gene_symbol": "C", "hpo_id": "HP:0001251"}])
        table(path / "edges.tsv", [{"subject": "HGNC:" + gene, "object": term,
                                    "predicate": "biolink:has_phenotype", "category": "biolink:GeneToPhenotypicFeatureAssociation"}
                                   for gene, term in [("A", "HP:0000648"), ("B", "HP:0007663"), ("C", "HP:0001251")]])
        table(path / "ann.tsv", [{"Gene_name": gene, "Annotation_mode": "split", "Exomiser_gene_pheno_score": score,
                                  "Human_pheno_evidence": gene + ":synthetic evidence"}
                                 for gene, score in [("A", "0.8"), ("B", "0.4"), ("C", "0")]])

    def compare(self, path, *extra, success=True):
        return self.command("compare", "--ranking", path / "rank.tsv", "--annotsv", path / "ann.tsv",
                            "--query", path / "query.tsv", "--phenotypes", path / "pheno.tsv", "--edges", path / "edges.tsv",
                            "--output", path / "comparison.tsv", "--summary-output", path / "summary.json", *extra, success=success)

    def test_shared_profile_keeps_only_core_ocular_terms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.prepare(path, "--enabled")
            terms = set((path / "query.txt").read_text().split(","))
            self.assertEqual(len(terms), 6)
            self.assertIn("HP:0000648", terms)
            self.assertNotIn("HP:0001251", terms)
            self.assertEqual({r["query_scope"] for r in read(path / "query.tsv")}, {"SHARED_HON_PROFILE"})

    def test_patient_file_overrides_shared_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            table(path / "patient.tsv", [{"HPO_ID": "HP:0001251"}])
            self.prepare(path, "--enabled", "--patient-hpo", path / "patient.tsv")
            self.assertEqual((path / "query.txt").read_text(), "HP:0001251")
            self.assertEqual(read(path / "query.tsv")[0]["query_scope"], "PATIENT_HPO_PROFILE")

    def test_invalid_ids_rejected_and_list_supported(self):
        self.assertEqual(hpo_ids(["HP:0000648", "HP:0000648"]), ["HP:0000648"])
        with self.assertRaises(ValueError):
            hpo_ids("optic atrophy")

    def test_backend_preflight_rejects_missing_or_incompatible_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            jar = path / "exomiser-cli-15.1.0.jar"
            with self.assertRaises(ValueError):
                check_backend_files(jar, path / "data", "2512", "2512")
            for resource in backend_resource_paths(jar, path / "data", "2512", "2512"):
                resource.parent.mkdir(parents=True, exist_ok=True)
                resource.write_text("synthetic fixture, not a real database")
            with self.assertRaisesRegex(ValueError, "valid ZIP/JAR"):
                check_backend_files(jar, path / "data", "2512", "2512")
            with zipfile.ZipFile(jar, "w") as handle:
                handle.writestr("META-INF/MANIFEST.MF", "Manifest-Version: 1.0\n")
            check_backend_files(jar, path / "data", "2512", "2512")
            with self.assertRaisesRegex(ValueError, "missing or empty"):
                check_backend_files(jar, path / "data", "2512", "2406")

    def native_fixture(self, path, java_major=21):
        jar = path / "exomiser-cli-15.1.0.jar"
        for resource in backend_resource_paths(jar, path / "data", "2512", "2512"):
            resource.parent.mkdir(parents=True, exist_ok=True)
            resource.write_text("synthetic interface fixture, not a real database")
        with zipfile.ZipFile(jar, "w") as handle:
            handle.writestr("META-INF/MANIFEST.MF", "Manifest-Version: 1.0\n")
        java = path / "fixture_java"
        # This executable tests our CLI interface and output parsing. It is
        # deliberately not described as an actual Exomiser scientific run.
        java.write_text("#!" + sys.executable + "\n" + '''import json, sys
from pathlib import Path
if '-version' in sys.argv:
    print('openjdk version "JAVA_MAJOR.0.1"', file=sys.stderr)
    raise SystemExit(0)
assert 'analyse' in sys.argv and '--analysis' in sys.argv
job = json.loads(Path(sys.argv[sys.argv.index('--analysis') + 1]).read_text())
assert job['preset'] == 'PHENOTYPE_ONLY'
assert job['sample']['genomeAssembly'] == 'hg38'
assert 'vcf' not in job['sample'] and 'ped' not in job['sample'] and 'sex' not in job['sample']
assert job['outputOptions']['numGenes'] == 0
assert job['outputOptions']['outputFormats'] == ['TSV_GENE']
out = Path(job['outputOptions']['outputDirectory']) / 'phenotype.genes.tsv'
out.write_text('#RANK\\tGENE_SYMBOL\\tEXOMISER_GENE_PHENO_SCORE\\tEXOMISER_GENE_COMBINED_SCORE\\tHUMAN_PHENO_EVIDENCE\\n'
               '1\\tA\\t0.8\\t0.01\\tA:synthetic evidence\\n'
               '2\\tB\\t0.4\\t0.99\\tB:synthetic evidence\\n'
               '3\\tC\\t0\\t0.9\\tC:synthetic evidence\\n'
               '4\\tNOT_IN_SV_RANKING\\t0.6\\t0.6\\tother\\n')
'''.replace("JAVA_MAJOR", str(java_major)))
        java.chmod(0o755)
        return jar, java

    def run_native(self, path, java_major=21, enabled=True, success=True):
        jar, java = self.native_fixture(path, java_major)
        native = path / "native"
        native.mkdir(exist_ok=True)
        arguments = ["run-offline", "--hpo-json", '["HP:0000648"]',
                     "--jar", jar, "--data-dir", path / "data", "--java", java,
                     "--genome-version", "2512", "--phenotype-version", "2512",
                     "--genes-output", native / "phenotype.genes.tsv",
                     "--job-output", native / "job.json", "--properties-output", native / "application.properties",
                     "--provenance-output", native / "run.json", "--log", native / "run.log"]
        if enabled:
            arguments.append("--enabled")
        return self.command(*arguments, success=success)

    def test_offline_cli_query_and_native_phenotype_comparison(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            before = (path / "rank.tsv").read_bytes()
            self.run_native(path)
            args = ["compare", "--ranking", path / "rank.tsv", "--exomiser-genes", path / "native/phenotype.genes.tsv",
                    "--backend-provenance", path / "native/run.json", "--query", path / "query.tsv",
                    "--phenotypes", path / "pheno.tsv", "--edges", path / "edges.tsv",
                    "--output", path / "comparison.tsv", "--summary-output", path / "summary.json", "--require-exomiser"]
            self.command(*args)
            rows = read(path / "comparison.tsv")
            self.assertEqual((path / "rank.tsv").read_bytes(), before)
            self.assertEqual(len(rows), 3)
            self.assertEqual([row["EXOMISER_GENE_PHENO_SCORE"] for row in rows], ["0.8", "0.4", "0.0"])
            self.assertEqual(rows[0]["EXOMISER_BACKEND"], "OFFLINE_CLI_PHENOTYPE_ONLY")
            self.assertEqual(rows[0]["EXOMISER_PHENOTYPE_DATA_VERSION"], "2512")
            self.assertEqual(rows[0]["EXOMISER_HUMAN_PHENO_EVIDENCE"], "A:synthetic evidence")
            properties = (path / "native/application.properties").read_text()
            self.assertIn("exomiser.hg38.data-version=2512", properties)
            self.assertNotIn("exomiser.hg19", properties)
            # A different saved query must not borrow scores from this run.
            self.prepare(path, "--enabled", "--shared-hpo-json", '"HP:0007663"')
            result = self.command(*args, success=False)
            self.assertIn("does not match the saved HPO query", result.stderr)
            self.prepare(path, "--enabled", "--shared-hpo-json", '"HP:0000648"')
            with (path / "native/phenotype.genes.tsv").open("a") as handle:
                handle.write("99\tNEW\t0.5\t0.5\tchanged\n")
            result = self.command(*args, success=False)
            self.assertIn("checksum differs", result.stderr)

    def test_offline_cli_rejects_java_8_and_disabled_needs_no_java(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            result = self.run_native(path, java_major=8, success=False)
            self.assertIn("requires Java 21", result.stderr)
            self.assertFalse((path / "native/phenotype.genes.tsv").exists())
            # Delete the Java executable; disabled mode must still succeed.
            (path / "fixture_java").unlink()
            self.command("run-offline", "--hpo-json", "[]", "--jar", "missing.jar", "--data-dir", "missing-data",
                         "--java", "missing-java", "--genome-version", "2512", "--phenotype-version", "2512",
                         "--genes-output", path / "native/phenotype.genes.tsv", "--job-output", path / "native/job.json",
                         "--properties-output", path / "native/application.properties", "--provenance-output", path / "native/run.json",
                         "--log", path / "native/run.log")
            self.assertEqual(json.loads((path / "native/run.json").read_text())["status"], "DISABLED")

    def test_full_and_multigene_rows_never_donate_scores(self):
        grouped = exomiser_by_gene([
            {"Gene_name": "A;B", "Annotation_mode": "full", "Exomiser_gene_pheno_score": "0.99"},
            {"Gene_name": "B", "Annotation_mode": "full", "Exomiser_gene_pheno_score": "0.99"},
            {"Gene_name": "A;B", "Annotation_mode": "split", "Exomiser_gene_pheno_score": "0.99"},
            {"Gene_name": "A", "Annotation_mode": "split", "Exomiser_gene_pheno_score": "0.2"},
        ])
        self.assertEqual(set(grouped), {"A"})
        self.assertEqual(extract_exomiser(grouped["A"], True, True)["EXOMISER_GENE_PHENO_SCORE"], 0.2)

    def test_invalid_negative_and_conflicting_scores_are_unavailable(self):
        for raw in ["-1", "NaN", "inf", "1.1", "."]:
            with self.subTest(raw=raw):
                result = extract_exomiser([{"Exomiser_gene_pheno_score": raw}], True, True)
                self.assertEqual(result["EXOMISER_GENE_PHENO_SCORE"], ".")
        result = extract_exomiser([{"Exomiser_gene_pheno_score": "0.8"}, {"Exomiser_gene_pheno_score": "0.2"}], True, True)
        self.assertEqual(result["EXOMISER_STATUS"], "CONFLICTING_GENE_SCORES")
        self.assertEqual(extract_exomiser([{"Exomiser_gene_pheno_score": "0"}], True, True)["EXOMISER_GENE_PHENO_SCORE"], 0)

    def test_comparison_preserves_original_ranking_and_retains_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            before = (path / "rank.tsv").read_bytes()
            self.compare(path, "--require-exomiser")
            self.assertEqual((path / "rank.tsv").read_bytes(), before)
            rows = read(path / "comparison.tsv")
            self.assertEqual([r["CUSTOM_GENE_DISPLAY_SCORE"] for r in rows], ["9.75", "6.25", "0.0"])
            self.assertEqual(rows[2]["EXOMISER_GENE_PHENO_SCORE"], "0.0")
            self.assertEqual(rows[0]["EXOMISER_HUMAN_PHENO_EVIDENCE"], "A:synthetic evidence")
            summary = json.loads((path / "summary.json").read_text())
            self.assertEqual(summary["paired_genes"], 3)
            self.assertEqual(summary["exomiser_scored_genes"], 3)
            self.assertEqual(summary["source_files"]["original_ranking"]["sha256"].__len__(), 64)

    def test_enabled_without_scores_fails_instead_of_silently_disabling(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            table(path / "ann.tsv", [{"Gene_name": "A", "Annotation_mode": "split", "Exomiser_gene_pheno_score": "-1"}])
            result = self.compare(path, "--require-exomiser", success=False)
            self.assertIn("no valid gene-specific scores", result.stderr)
            self.assertFalse((path / "comparison.tsv").exists())

    def test_missing_gene_annotations_do_not_become_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            table(path / "pheno.tsv", [{"gene_symbol": "A", "hpo_id": "HP:0000648"}])
            self.compare(path, "--require-exomiser")
            self.assertEqual(read(path / "comparison.tsv")[1]["CUSTOM_QUERY_PHENOTYPE_SCORE_0_1"], ".")
            self.assertEqual(json.loads((path / "summary.json").read_text())["paired_genes"], 1)

    def test_query_missing_from_monarch_prevents_false_comparison(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            self.prepare(path, "--enabled", "--shared-hpo-json", '"HP:9999999"')
            self.compare(path)
            self.assertEqual(read(path / "comparison.tsv")[0]["CUSTOM_QUERY_STATUS"], "QUERY_TERMS_MISSING_FROM_MONARCH")
            self.assertEqual(json.loads((path / "summary.json").read_text())["paired_genes"], 0)

    def test_disabled_and_empty_comparisons_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            self.prepare(path)
            self.compare(path)
            self.assertEqual(read(path / "comparison.tsv")[0]["EXOMISER_STATUS"], "DISABLED")
            table(path / "rank.tsv", [], ["gene", "GENE_RANK_GLOBAL"])
            self.compare(path)
            self.assertEqual(read(path / "comparison.tsv"), [])

    def test_plot_exports_readable_comparison_and_source_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            self.fixtures(path)
            self.compare(path)
            result = subprocess.run([sys.executable, str(ROOT / "plots/plot_gene_hpo_heatmap.py"),
                                     "--comparison-only", "--exomiser-comparison", str(path / "comparison.tsv"),
                                     "--out-prefix", str(path / "figure")], capture_output=True, text=True,
                                    env={**os.environ, "MPLBACKEND": "Agg"})
            self.assertEqual(result.returncode, 0, result.stderr)
            for suffix in [".png", ".pdf", ".svg", "_source.tsv"]:
                self.assertTrue((path / ("figure" + suffix)).is_file())


if __name__ == "__main__":
    unittest.main()
