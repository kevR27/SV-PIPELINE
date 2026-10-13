"""Focused checks for failures/reuse found in the October 2026 input audit."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import check_run_cache as cache
import build_snv_sv_candidates as pairs
import annotate_alphagenome_atlas as atlas

spec = importlib.util.spec_from_file_location("srs_samplot", ROOT / "scripts/run_samplot_candidates.py")
samplot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(samplot)


class PipelineIOAuditTests(unittest.TestCase):
    def fixture(self, path):
        backend = path / "backend with spaces"
        backend.mkdir()
        (backend / "reference.tsv").write_text("reference version one\n")
        query = path / "query.vcf.gz"
        with gzip.open(query, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        output = path / "result.vcf.gz"
        with gzip.open(output, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        current = cache.snapshot([backend], [query], {"mode": "annotate --all"})
        manifest = path / "completed.json"
        manifest.write_text(json.dumps({"snapshot": current, "outputs": cache.output_records([output])}))
        return backend, query, output, current, manifest

    def test_unsigned_or_changed_backend_outputs_are_not_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            backend, query, output, current, manifest = self.fixture(path)
            self.assertTrue(cache.matches(manifest, current))
            self.assertFalse(cache.matches(path / "missing.json", current))
            (backend / "reference.tsv").write_text("a different population database\n")
            changed = cache.snapshot([backend], [query], {"mode": "annotate --all"})
            self.assertFalse(cache.matches(manifest, changed))

    def test_corrupt_or_changed_generated_outputs_are_not_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            _, _, output, current, manifest = self.fixture(Path(temp))
            output.write_bytes(output.read_bytes()[:-8])
            self.assertFalse(cache.matches(manifest, current))
            with gzip.open(output, "wt") as handle:
                handle.write("different result\n")
            self.assertFalse(cache.matches(manifest, current))

    def test_needlr_unquoted_internal_paths_are_rejected_early(self):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/check_run_cache.py"),
                                 "snapshot", "--tool-paths", "/example/path with spaces"], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("needLR 4.1 uses unquoted paths", result.stderr)

    def test_same_uncompressed_query_and_changed_settings(self):
        with tempfile.TemporaryDirectory() as temp:
            backend, query, _, current, manifest = self.fixture(Path(temp))
            with gzip.open(query, "rt") as handle:
                text = handle.read()
            with gzip.open(query, "wt", compresslevel=1) as handle:
                handle.write(text)
            self.assertTrue(cache.matches(manifest, cache.snapshot([backend], [query], {"mode": "annotate --all"})))
            self.assertFalse(cache.matches(manifest, cache.snapshot([backend], [query], {"mode": "different"})))

    def test_empty_straglr_is_a_valid_negative_but_bad_bed_is_not(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            columns = "chrom start end target_repeat locus coverage genotype actual_repeat read_name copy_number size read_status".split()
            (path / "repeat.tsv").write_text("#" + "\t".join(columns) + "\n")
            (path / "repeat.bed").write_text("")
            command = [sys.executable, str(ROOT / "scripts/annotate_straglr.py"), "--straglr-tsv", str(path / "repeat.tsv"),
                       "--straglr-bed", str(path / "repeat.bed"), "--gene-bed", str(path / "unused.bed"), "--output", str(path / "out.tsv")]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("NO_REPEAT_LOCI", result.stdout)
            with (path / "out.tsv").open() as handle:
                self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), [])
            (path / "repeat.bed").write_text("chr1\tinvalid\t100\n")
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)

    def test_srs_samplot_cannot_count_a_stale_image_as_success(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            (path / "in.tsv").write_text("SV_ID\tCHROM\tSTART\tEND\tSVTYPE\tGENE\nv\tchr1\t100\t200\tDEL\tOPA1\n")
            (path / "images").mkdir()
            stale = path / "images/001_DEL_chr1_100_200_v_OPA1.png"
            stale.write_bytes(b"old image")
            argv = ["samplot", "--input", str(path / "in.tsv"), "--sample", "SYNTHETIC", "--bam", "not-real.bam",
                    "--reference", "not-real.fa", "--out-dir", str(path / "images"), "--manifest", str(path / "manifest.tsv")]
            with patch.object(sys, "argv", argv), patch.object(samplot.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout="")):
                self.assertEqual(samplot.main(), 1)
            self.assertFalse(stale.exists())
            with (path / "manifest.tsv").open() as handle:
                self.assertEqual(next(csv.DictReader(handle, delimiter="\t"))["STATUS"], "FAILED")

    def test_genomewide_vep_master_is_streamed_into_existing_focused_view(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "full.tsv.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("#Uploaded_variation\tLocation\tAllele\tConsequence\tExtra\n"
                             "1_100_A/G\t1:100\tG\tintron_variant\tSYMBOL=OPA1;IMPACT=MODIFIER\n"
                             "1_101_A/G\t1:101\tG\tmissense_variant\tSYMBOL=OTHER;IMPACT=MODERATE\n"
                             "1_102_A/G\t1:102\tG\tsplice_donor_variant\tSYMBOL=OPA1;IMPACT=HIGH\n")
            self.assertEqual(len(pairs.read_vep(str(path))), 3)
            focused = pairs.read_vep(str(path), {"OPA1"})
            self.assertEqual(focused["Uploaded_variation"].tolist(), ["1_102_A/G"])

    def test_atlas_matches_single_base_vep_range(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "context.tsv"
            path.write_text("#Uploaded_variation\tLocation\tAllele\tConsequence\tExtra\n"
                            "v\t1:100-100\tG\tintron_variant\tUPLOADED_ALLELE=A/G;SYMBOL=OPA1\n")
            row = {"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "G"}
            atlas.add_vep_context([row], path)
            self.assertEqual(row["VEP_CONTEXT_STATUS"], "EXACT_SNV_ALLELE_MATCH")

    def test_runner_and_source_output_contracts(self):
        lrs = (ROOT / "snakemake_pipelines/lrs/Snakefile_LRS_update").read_text()
        self.assertNotIn("cp {params.outbase}.table.txt {output.table}", lrs)
        self.assertNotIn("Bootstrap existing pre-checksum", lrs)
        srs = (ROOT / "snakemake_pipelines/srs_wgs_pipeline/Snakefile_SRS_WGS").read_text()
        self.assertIn("/opt/deepvariant/bin/run_deepvariant", srs)
        self.assertIn("apptainer exec --bind", srs)
        self.assertIn("--cache_version {VEP_CACHE_VERSION}", srs)


if __name__ == "__main__":
    unittest.main()
