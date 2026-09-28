#!/usr/bin/env python3
"""Reference conversion and installation tests using synthetic transcripts."""
import argparse
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import rebuild_annotsv_gene_resource as rebuild
import install_annotsv_gene_resource as installer


def fixtures(folder, symbol="", coding=True):
    tx, gene = "ENST00000000001.2", "ENSG00000000001.3"
    gp = folder / "source.genePred"
    info = folder / "info.tsv"
    gp_fields = [tx, "1", "+", "100", "200", "120" if coding else "200",
                 "180" if coding else "200", "1", "100,", "200,", "0", gene,
                 "cmpl" if coding else "none", "cmpl" if coding else "none", "0,"]
    gp.write_text("\t".join(gp_fields) + "\n")
    info_fields = [tx, gene, "ensembl", "1", "100", "200", "+", "", symbol,
                   "", "protein_coding" if coding else "lncRNA", "protein_coding" if coding else "lncRNA"]
    info.write_text("#" + "\t".join(rebuild.INFO_FIELDS) + "\n" + "\t".join(info_fields) + "\n")
    return gp, info


class GeneRebuildTests(unittest.TestCase):
    def test_empty_name_retains_gene_id_and_transcript_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, identities, skipped = rebuild.convert_tables(*fixtures(Path(tmp)))
            tx = "ENST00000000001"
            self.assertEqual(rows[tx][4], "ENSG00000000001")
            self.assertEqual(identities[tx][1], "2")
            self.assertEqual(identities[tx][3], "3")
            self.assertEqual(identities[tx][6], "ENSEMBL_GENE_ID_NO_SYMBOL")
            self.assertEqual(skipped, {})

    def test_unk_is_a_valid_gene_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, _, _ = rebuild.convert_tables(*fixtures(Path(tmp), symbol="UNK"))
            self.assertEqual(next(iter(rows.values()))[4], "UNK")

    def test_unnamed_noncoding_transcript_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, _, _ = rebuild.convert_tables(*fixtures(Path(tmp), coding=False))
            row = next(iter(rows.values()))
            self.assertEqual(row[4], "ENSG00000000001")
            self.assertEqual(row[6], row[7])

    def test_mismatched_gene_metadata_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            gp, info = fixtures(Path(tmp))
            info.write_text(info.read_text().replace("ENSG00000000001.3", "ENSG00000000002.1"))
            with self.assertRaisesRegex(ValueError, "Conflicting source gene"):
                rebuild.convert_tables(gp, info)

    def test_duplicate_or_missing_transcript_information_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            gp, info = fixtures(Path(tmp))
            original = gp.read_text()
            gp.write_text(original + original)
            with self.assertRaisesRegex(ValueError, "Duplicate genePred"):
                rebuild.convert_tables(gp, info)
            gp.write_text(original.replace("ENST00000000001", "ENST00000000002"))
            with self.assertRaisesRegex(ValueError, "Missing source information"):
                rebuild.convert_tables(gp, info)

    def test_mitochondrial_names_normalize_and_other_contigs_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            gp, info = fixtures(Path(tmp))
            fields = gp.read_text().rstrip("\n").split("\t")
            fields[1] = "MT"
            gp.write_text("\t".join(fields) + "\n")
            header, entry = info.read_text().splitlines()
            fields = entry.split("\t")
            fields[3] = "MT"
            info.write_text(header + "\n" + "\t".join(fields) + "\n")
            rows, _, _ = rebuild.convert_tables(gp, info)
            self.assertEqual(next(iter(rows.values()))[0], "M")
            g2 = gp.read_text().replace("ENST00000000001", "ENST00000000002").replace("\tMT\t", "\tGL000001.1\t")
            i2 = info.read_text().splitlines()[1].replace("ENST00000000001", "ENST00000000002").replace("\tMT\t", "\tGL000001.1\t")
            gp.write_text(gp.read_text() + g2)
            info.write_text(info.read_text() + i2 + "\n")
            rows, _, skipped = rebuild.convert_tables(gp, info)
            self.assertEqual(len(rows), 1)
            self.assertEqual(skipped, {"GL000001.1": 1})

    def test_comparison_distinguishes_name_repair_from_coordinate_changes(self):
        old = {"ENST00000000001": ["1", "100", "200", "+", "cmpl", "ENST00000000001", "120", "180", "100,", "200,"]}
        new = {tx: row.copy() for tx, row in old.items()}
        new["ENST00000000001"][4] = "ENSG00000000001"
        counts, _ = rebuild.compare(old, new)
        self.assertEqual(counts, {"INVALID_NAME_REPLACED": 1})
        new["ENST00000000001"][6] = "121"
        counts, _ = rebuild.compare(old, new)
        self.assertEqual(counts["TRANSCRIPT_GEOMETRY_CHANGED"], 1)

    def test_wrong_genome_build_is_rejected_for_gzip_and_plain_gtf(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            for suffix, opener in ((".gtf", open), (".gtf.gz", gzip.open)):
                p = d / ("source" + suffix)
                with opener(p, "wt") as f:
                    f.write("#!genome-build GRCh37\n")
                with self.assertRaisesRegex(ValueError, "GRCh38"):
                    rebuild.copy_gtf(p, d / "unpacked.gtf")

    def make_build(self, folder, change_geometry=False):
        gp, info = fixtures(folder)
        old = folder / "old.bed"
        old.write_text("1\t100\t200\t+\tcmpl\tENST00000000001\t120\t180\t100,\t200,\n")
        if change_geometry:
            gp.write_text(gp.read_text().replace("\t120\t", "\t121\t"))
        gtf = folder / "source.gtf"
        gtf.write_text("#!genome-build GRCh38.p14\n")
        args = argparse.Namespace(gtf=gtf, existing_bed=old, outdir=folder / "build",
                                  gtf_to_genepred=sys.executable, source_release=113,
                                  allow_reference_update=False)
        def conversion(command, **kwargs):
            Path(command[-1]).write_text(gp.read_text())
            Path(next(v.split("=", 1)[1] for v in command if v.startswith("-infoOut="))).write_text(info.read_text())
        with patch.object(rebuild.subprocess, "run", side_effect=conversion):
            rebuild.rebuild(args)
        return args.outdir, old

    def test_incompatible_source_stops_before_writing_installable_bed(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            with self.assertRaisesRegex(ValueError, "differs beyond"):
                self.make_build(d, change_geometry=True)
            self.assertFalse((d / "build" / installer.BED).exists())
            report = json.loads((d / "build/reference_build.json").read_text())
            self.assertEqual(report["status"], "FAILED")

    def test_install_backs_up_originals_and_refreshes_only_ensembl_promoters(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            build, old = self.make_build(d)
            root = d / "installed"
            (root / installer.GENES).mkdir(parents=True)
            (root / installer.BED).write_bytes(old.read_bytes())
            (root / installer.VERSIONS).write_text("ENST00000000001\t1\n")
            (root / installer.REGULATORY).mkdir(parents=True)
            promoter = root / installer.REGULATORY / "promoter_500bp_ENSEMBL_GRCh38.sorted.bed"
            other = root / installer.REGULATORY / "promoter_500bp_RefSeq_GRCh38.sorted.bed"
            promoter.write_text("old promoter\n")
            other.write_text("keep this\n")
            (root / ".installed").write_text("")
            backup = installer.install(build, root)
            self.assertEqual((backup / installer.BED).read_bytes(), old.read_bytes())
            self.assertEqual((backup / installer.VERSIONS).read_text(), "ENST00000000001\t1\n")
            self.assertFalse(promoter.exists())
            self.assertFalse((root / ".installed").exists())
            self.assertEqual(other.read_text(), "keep this\n")
            self.assertIn("ENSG00000000001", (root / installer.BED).read_text())
            self.assertEqual((root / installer.VERSIONS).read_text(), "ENST00000000001\t2\n")

    def test_install_refuses_changed_source_and_changed_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            build, old = self.make_build(d)
            root = d / "installed"
            (root / installer.GENES).mkdir(parents=True)
            (root / installer.BED).write_text(old.read_text().replace("cmpl", "OPA1"))
            with self.assertRaisesRegex(ValueError, "installed BED differs"):
                installer.install(build, root)
            (root / installer.BED).write_bytes(old.read_bytes())
            (build / installer.BED).write_text((build / installer.BED).read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "changed after validation"):
                installer.install(build, root)

    def test_install_restores_original_if_second_file_install_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            build, old = self.make_build(d)
            root = d / "installed"
            (root / installer.GENES).mkdir(parents=True)
            (root / installer.BED).write_bytes(old.read_bytes())
            (root / installer.VERSIONS).write_text("ENST00000000001\t1\n")
            real_replace = installer.os.replace
            count = 0
            def fail_second(source, destination):
                nonlocal count
                count += 1
                if count == 2:
                    raise OSError("simulated write failure")
                return real_replace(source, destination)
            with patch.object(installer.os, "replace", side_effect=fail_second):
                with self.assertRaisesRegex(OSError, "simulated write failure"):
                    installer.install(build, root)
            self.assertEqual((root / installer.BED).read_bytes(), old.read_bytes())
            self.assertEqual((root / installer.VERSIONS).read_text(), "ENST00000000001\t1\n")
            self.assertFalse((root / ".gene-reference-update.lock").exists())


if __name__ == "__main__":
    unittest.main()
