"""Focused allele/coordinate checks; synthetic data, no patient/API analysis."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import annotate_alphagenome_atlas as atlas


def table(path):
    with atlas.open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AtlasTests(unittest.TestCase):
    def fixture(self, directory):
        path = Path(directory)
        (path / "ref.fa").write_text(">chr1\n" + "A" * 2000000 + "\n")
        pysam.faidx(str(path / "ref.fa"))
        header = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=2000000>\n"
                  '##FILTER=<ID=LowQual,Description="Low quality">\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                  '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n')
        (path / "small.vcf").write_text(header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n"
            "chr1\t10000\t.\tA\tC,G\t50\tPASS\t.\tGT:GQ:DP\t0/2:40:30\n"
            "chr1\t10001\t.\tA\tT\t50\tLowQual\t.\tGT:GQ:DP\t0/1:10:4\n"
            "chr1\t10002\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30\n"
            "chr1\t10003\t.\tC\tG\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30\n"
            "chr1\t10004\t.\tA\tAT\t50\tPASS\t.\tGT:GQ:DP\t0/1:40:30\n")
        svheader = (header + '##INFO=<ID=END,Number=1,Type=Integer,Description="END">\n'
                    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">\n'
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">\n'
                    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="IMPRECISE">\n'
                    '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CIPOS">\n')
        (path / "sv.vcf").write_text(svheader + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n"
            "chr1\t600000\tdel\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=650000;SVLEN=-50000\tGT:GQ:DP\t0/1:40:30\n"
            "chr1\t700000\tdup\tA\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=750000;SVLEN=50000\tGT:GQ:DP\t0/1:40:30\n"
            "chr1\t800000\timprecise\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=850000;SVLEN=-50000;IMPRECISE\tGT:GQ:DP\t0/1:40:30\n"
            "chr1\t900000\tbad_len\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=950000;SVLEN=-50001\tGT:GQ:DP\t0/1:40:30\n")
        (path / "avi.tsv").write_text("#chrom\tpos\tref\talt\traw_score\tphred\n"
            "chr1\t10000\tA\tC\t-1.2\t0\nchr1\t10000\tA\tG\t1.5\t25\n"
            "chr1\t10001\tA\tT\t1.8\t30\n")
        pysam.tabix_compress(str(path / "avi.tsv"), str(path / "avi.tsv.gz"), force=True)
        pysam.tabix_index(str(path / "avi.tsv.gz"), seq_col=0, start_col=1, end_col=1, force=True)
        (path / "genes.bed").write_text("chr1\t9900\t11000\tOPA1\n")
        (path / "panel.txt").write_text("OPA1\n")
        (path / "vep.tsv").write_text("#Uploaded_variation\tLocation\tAllele\tGene\tConsequence\tExtra\n"
            "chr1_10000_A/C,G\t1:10000\tG\tENSG1\tintron_variant\tSYMBOL=OPA1;UPLOADED_ALLELE=A/C,G\n"
            "chr1_10001_A/T\t1:10001\tT\tENSG2\tregulatory_region_variant\tUPLOADED_ALLELE=A/T\n")
        return path

    def arguments(self, path):
        return SimpleNamespace(vcf=str(path / "small.vcf"), sv_vcf=str(path / "sv.vcf"), reference=str(path / "ref.fa"),
            avi=str(path / "avi.tsv.gz"), gene_bed=str(path / "genes.bed"), panel=str(path / "panel.txt"), vep=str(path / "vep.tsv"),
            sample="SYNTHETIC", release="SYNTHETIC_FIXTURE", assembly="GRCh38", dna_tissue="muscle", bin_size=1000,
            min_gq=20, min_dp=10, review_phred=20, output=str(path / "alleles.tsv.gz"), review_output=str(path / "review.tsv"),
            sv_output=str(path / "scope.tsv"), summary=str(path / "summary.json"))

    def test_all_alleles_retained_and_missing_not_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(directory)
            atlas.annotate(self.arguments(path))
            rows = table(path / "alleles.tsv.gz")
            self.assertEqual(len(rows), 6)
            self.assertEqual((rows[0]["AVI_RAW"], rows[0]["AVI_PHRED"]), ("-1.2", "0.0"))
            self.assertEqual(rows[0]["ALT_PRESENT"], "NO")
            self.assertEqual(rows[1]["ALT_PRESENT"], "YES")
            self.assertIn("FILTER_NOT_PASS", rows[2]["QC_FLAGS"])
            self.assertEqual(rows[3]["AVI_STATUS"], "ALLELE_NOT_IN_ATLAS")
            self.assertEqual(rows[3]["AVI_PHRED"], ".")
            self.assertEqual(rows[4]["AVI_STATUS"], "REF_MISMATCH")
            self.assertEqual(rows[5]["AVI_STATUS"], "NOT_SUPPORTED_BY_OFFLINE_SNV_ATLAS")
            review = table(path / "review.tsv")
            self.assertEqual(len(review), 2)
            self.assertIn("regulatory_region_variant", review[0]["VEP_CONSEQUENCES"])
            self.assertIn("intron_variant", review[1]["VEP_CONSEQUENCES"])
            self.assertTrue(all(r["AVI_STATUS"] == "SV_NOT_SUPPORTED_BY_OFFLINE_SNV_ATLAS" for r in table(path / "scope.tsv")))
            data = json.loads((path / "summary.json").read_text())
            self.assertEqual(sum(data["alleles_by_status"].values()), 6)
            self.assertEqual(data["dna_sample_tissue"], "muscle")

    def test_bad_resource_and_conflicts_do_not_silently_score(self):
        for record in ["chr1\t3\tA\tG\tNaN\t20", "chr1\t3\tA\tG\t1\t-1", "chr1\t3\tA\tG\t1\t2\tother", "chr1\t3\tAA\tG\t1\t20"]:
            with self.subTest(record=record), self.assertRaises(ValueError):
                atlas.parse_avi(record)
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(directory)
            (path / "bad.tsv").write_text("chr1\t10000\tA\tG\t1\t20\nchr1\t10000\tA\tG\t2\t20\n")
            pysam.tabix_compress(str(path / "bad.tsv"), str(path / "bad.tsv.gz"), force=True)
            pysam.tabix_index(str(path / "bad.tsv.gz"), seq_col=0, start_col=1, end_col=1, force=True)
            lookup = atlas.AviLookup(path / "bad.tsv.gz")
            with self.assertRaisesRegex(ValueError, "Conflicting"):
                lookup.get("1", 10000, "A", "G")
            lookup.tabix.close()

    def test_vep_different_ref_is_not_a_gene_only_match(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(directory)
            (path / "vep.tsv").write_text("#Uploaded_variation\tLocation\tAllele\tConsequence\tExtra\n"
                "wrong\t1:10000\tG\tmissense_variant\tUPLOADED_ALLELE=C/G;SYMBOL=WRONG\n")
            row = {"CHROM": "chr1", "POS": 10000, "REF": "A", "ALT": "G"}
            atlas.add_vep_context([row], path / "vep.tsv")
            self.assertEqual(row["VEP_CONTEXT_STATUS"], "EXACT_ALLELE_NOT_FOUND")

    def test_zero_based_index_is_rejected_and_caller_columns_are_not_patient_gt(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(directory)
            pysam.tabix_index(str(path / "avi.tsv.gz"), seq_col=0, start_col=1, end_col=1, zerobased=True, force=True)
            with self.assertRaisesRegex(ValueError, "1-based"):
                atlas.AviLookup(path / "avi.tsv.gz")
            source = (path / "sv.vcf").read_text()
            source = source.replace("FORMAT\tSYNTHETIC", "FORMAT\tcaller_one\tcaller_two")
            source = source.replace("GT:GQ:DP\t0/1:40:30", "GT:GQ:DP\t0/1:40:30\t./.:.:.")
            (path / "sv.vcf").write_text(source)
            rows = list(atlas.sv_inventory(path / "sv.vcf", "SYNTHETIC"))
            self.assertEqual(rows[0]["GT"], ".")
            self.assertEqual(rows[0]["GENOTYPE_SCOPE"], "MERGED_CALLER_COLUMNS_UNRESOLVED")
            self.assertIn("END_SVLEN_CONFLICT", rows[-1]["BREAKPOINT_FLAGS"])

    def test_large_del_anchor_coordinates_and_required_review(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(directory)
            with pysam.VariantFile(str(path / "sv.vcf")) as vcf, pysam.FastaFile(str(path / "ref.fa")) as fasta:
                records = list(vcf)
                selection = {"SOURCE": "SV", "CONFIRMED_BREAKPOINTS": "YES", "CONFIRMED_ALT_PRESENT": "YES"}
                chrom, pos, ref, alt, start, end, label = atlas.selected_allele(records[0], 1, selection, fasta, 1048576, 10000)
                self.assertEqual((chrom, pos, len(ref), alt), ("chr1", 600000, 50001, "A"))
                self.assertEqual(len(ref) - len(alt), 50000)
                self.assertLess(start, 599999)
                self.assertGreater(end, 650000)
                self.assertEqual(label, "EXPLORATORY_DEL_UNVALIDATED_NO_AVI")
                for rec in records[1:]:
                    with self.assertRaises(ValueError):
                        atlas.selected_allele(rec, 1, selection, fasta, 1048576, 10000)
                with self.assertRaisesRegex(ValueError, "manual review"):
                    atlas.selected_allele(records[0], 1, {"SOURCE": "SV"}, fasta, 1048576, 10000)
                with self.assertRaisesRegex(ValueError, "context"):
                    atlas.selected_allele(records[0], 1, selection, fasta, 16384, 1000)

    def test_overview_and_native_sdk_plots_use_cached_source_tracks(self):
        from alphagenome.data import genome, track_data
        from alphagenome.models import dna_client
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(directory)
            atlas.annotate(self.arguments(path))
            atlas.summary_plot(SimpleNamespace(summary=str(path / "summary.json"), outdir=str(path / "plots")))
            self.assertTrue((path / "plots/atlas_overview.source.tsv").is_file())
            (path / "selection.tsv").write_text("SOURCE\tCHROM\tPOS\tREF\tALT\tSV_ID\nSMALL\tchr1\t10000\tA\tG\t.\n")
            args = SimpleNamespace(vcf=str(path / "small.vcf"), sv_vcf=str(path / "sv.vcf"), reference=str(path / "ref.fa"),
                sample="SYNTHETIC", selection=str(path / "selection.tsv"), outdir=str(path / "plots/native"),
                cache_dir=str(path / "cache"), manifest=str(path / "plots/native_manifest.tsv"),
                ontology_terms_json='["UBERON:0000178"]', max_variants=10, min_gq=20, min_dp=10,
                sequence_length=16384, min_flank=1000, model_version="ALL_FOLDS", cache_only=False)
            def prediction(interval, variant, **kwargs):
                metadata = pd.DataFrame({"name": ["SYNTHETIC blood RNA"], "strand": ["+"], "ontology_curie": ["UBERON:0000178"]})
                ref = np.ones((interval.width // 128, 1), dtype=np.float32)
                alt = ref.copy()
                alt[len(alt) // 2:] = 0.5
                return SimpleNamespace(reference=SimpleNamespace(rna_seq=track_data.TrackData(ref, metadata, resolution=128, interval=interval), splice_sites=None),
                                       alternate=SimpleNamespace(rna_seq=track_data.TrackData(alt, metadata, resolution=128, interval=interval), splice_sites=None))
            client = SimpleNamespace(predict_variant=prediction)
            with patch.object(dna_client, "create", return_value=client) as create, patch.dict("os.environ", {"ALPHAGENOME_API_KEY": "SYNTHETIC_NOT_A_KEY"}):
                atlas.native_plots(args)
                create.assert_called_once()
            manifest = table(args.manifest)
            folder = Path(manifest[0]["PLOT_FOLDER"])
            for ext in ["pdf", "svg", "png"]:
                self.assertGreater((folder / f"rna_seq.{ext}").stat().st_size, 100)
            self.assertEqual(manifest[0]["LABEL"], "SMALL_VARIANT_RESEARCH")
            args.cache_only = True
            with patch.object(dna_client, "create", side_effect=AssertionError("Cache should avoid API")):
                atlas.native_plots(args)
            args.ontology_terms_json = '["UBERON:0000955"]'
            with self.assertRaisesRegex(ValueError, "No cached tracks"):
                atlas.native_plots(args)

    def test_both_workflows_keep_original_callers_and_shared_module(self):
        for folder, snake, caller, context in [
            ("lrs", "Snakefile_LRS_postprocess", "snp_clair3/{sample}.vcf.gz", "longphase.vep.txt"),
            ("srs_wgs_pipeline", "Snakefile_SRS_WGS", "snp_deepvariant/{sample}.vcf.gz", "genomewide.vep.tsv.gz")]:
            text = (ROOT / "snakemake_pipelines" / folder / snake).read_text()
            self.assertIn(caller, text)
            self.assertIn(context, text)
            self.assertIn('include: "../common/alphagenome_atlas.smk"', text)
            self.assertIn("atlas=lambda wc: ATLAS_ANNOTATION_TARGETS", text)


if __name__ == "__main__":
    unittest.main()
