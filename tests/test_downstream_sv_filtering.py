"""Boundary and evidence-semantics tests; all variants are synthetic."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "plots"))
import filter_sv_candidates_for_review as filtering
import plot_downstream_sv_filtering as plotting


def base_row():
    return {"SV_ID": "synthetic1", "GENE": "TEST_GENE", "CHROM": "chr1",
            "START": "10000", "END": "12000", "SVTYPE": "DEL", "SVLEN": "-2000",
            "FINAL_EVENT_RANK_GLOBAL": "7", "PANEL_STATUS": "PANEL_GENE",
            "CALLER_COUNT": "3", "CALLERS": "Sniffles2;cuteSV;Delly",
            "CALLER_READ_SUPPORT": "Sniffles2:3;cuteSV:4;Delly:5",
            "NEEDLR_AF": "0.01", "VEP_MATCH_STATUS": "GENE_TRANSCRIPT_MATCH",
            "VEP_TRANSCRIPT_REGION_CLASS": "CODING"}


class DownstreamFilteringTests(unittest.TestCase):
    def annotate(self, **updates):
        row = {**base_row(), **updates}
        return filtering.annotate_row(row, "SYNTHETIC", 1, 1000, 0.01, 3, 5)

    def run_filter(self, source, out, *extra):
        with patch.object(sys, "argv", ["filter", "--input", str(source), "--out-dir", str(out), "--sample", "SYNTHETIC", *extra]), contextlib.redirect_stdout(io.StringIO()):
            filtering.main()
        return pd.read_csv(out / "all_records.annotated.tsv", sep="\t", dtype=str, keep_default_na=False)

    def test_size_cutoff_and_signed_lengths(self):
        for size, expected in [(999, "LT_CUTOFF"), (1000, "AT_CUTOFF"), (1001, "GT_CUTOFF"), (10_000_000, "GT_CUTOFF")]:
            with self.subTest(size=size):
                row = self.annotate(SVLEN=str(-size))
                self.assertEqual(row["DOWNSTREAM_SIZE_GROUP"], expected)
                self.assertEqual(row["DOWNSTREAM_SIZE_BP"], size)

    def test_insertion_and_breakend_lengths_are_not_invented(self):
        row = {**base_row(), "SVTYPE": "INS", "SVLEN": ".", "SV_SPAN_BP": "5000", "END": "50000"}
        self.assertIsNone(filtering.event_size(row)[0])
        row.update(REF="N", ALT="N" + "A"*1001)
        self.assertEqual(filtering.event_size(row), (1001, "LITERAL_INSERTED_SEQUENCE"))
        self.assertIsNone(filtering.event_size({**base_row(), "SVTYPE": "BND", "SVLEN": "100000"})[0])

    def test_y_and_mt_are_removed_from_either_endpoint(self):
        for update in [{"CHROM": "chrY"}, {"CHR2": "Y"}, {"ALT": "N[chrY:100["}]:
            self.assertIn("EXCLUDE_Y", self.annotate(**update)["DOWNSTREAM_EXCLUSION_REASONS"])
        for update in [{"CHROM": "chrM"}, {"CHROM": "MT"}, {"ALT": "N]chrMT:100]"}]:
            self.assertIn("EXCLUDE_MITOGENOME", self.annotate(**update)["DOWNSTREAM_EXCLUSION_REASONS"])
        self.assertEqual(self.annotate(MITOCARTA="NUCLEAR_ENCODED")["DOWNSTREAM_CHROMOSOME_SCOPE"], "NUCLEAR_WITHOUT_Y")

    def test_af_boundary_zero_and_unknown_are_distinct(self):
        for af in ["0", "0.001", "0.01"]:
            row = self.annotate(NEEDLR_AF=af)
            self.assertEqual(row["DOWNSTREAM_AF_STATUS"], "MEASURED_AT_OR_BELOW_CUTOFF")
            self.assertEqual(row["DOWNSTREAM_EXCLUSION_REASONS"], ".")
        self.assertIn("AF_ABOVE_CUTOFF", self.annotate(NEEDLR_AF="0.01001")["DOWNSTREAM_EXCLUSION_REASONS"])
        row = self.annotate(NEEDLR_AF=".", AF="0.5", ALLELE_AF="0.5")
        self.assertEqual(row["DOWNSTREAM_AF_STATUS"], "AF_UNKNOWN")
        self.assertEqual(row["DOWNSTREAM_MAX_AF"], ".")

    def test_highest_population_af_wins_and_invalid_values_are_reviewed(self):
        row = self.annotate(NEEDLR_AF="0.001", GNOMAD_SV_AF="0.004;0.2", GNOMAD_SV_EXACT_MATCH="YES")
        self.assertEqual(row["DOWNSTREAM_MAX_AF"], 0.2)
        self.assertIn("AF_ABOVE_CUTOFF", row["DOWNSTREAM_EXCLUSION_REASONS"])
        for af in ["1.5", "-0.1", "garbage", "0.001;garbage"]:
            self.assertEqual(self.annotate(NEEDLR_AF=af)["DOWNSTREAM_AF_STATUS"], "INVALID_AF_REVIEW")
        row = self.annotate(NEEDLR_AF="0.8", NEEDLR_STATUS="NO_MATCH", GNOMAD_SV_AF="0.4", GNOMAD_SV_EXACT_MATCH="NO")
        self.assertEqual(row["DOWNSTREAM_AF_STATUS"], "AF_UNKNOWN")

    def test_actual_vep_region_labels_and_mixed_transcripts(self):
        for fields, expected in [({"VEP_TRANSCRIPT_REGION_CLASS": "CODING"}, "EXONIC_OR_SPLICE"),
                                 ({"VEP_TRANSCRIPT_REGION_CLASS": "INTRONIC"}, "INTRONIC"),
                                 ({"VEP_TRANSCRIPT_REGION_CLASS": "CODING;INTRONIC"}, "EXONIC_AND_INTRONIC"),
                                 ({"VEP_TRANSCRIPT_REGION_CLASS": ".", "VEP_CONSEQUENCES": "splice_donor_variant"}, "EXONIC_OR_SPLICE"),
                                 ({"VEP_TRANSCRIPT_REGION_CLASS": "TRANSCRIPT_ABLATION"}, "EXONIC_OR_SPLICE")]:
            self.assertEqual(self.annotate(**fields)["DOWNSTREAM_CONTEXT"], expected)
        self.assertEqual(self.annotate(VEP_TRANSCRIPT_REGION_CLASS=".", SV_GENE_EFFECT="PARTIAL_GENE_DELETION")["DOWNSTREAM_CONTEXT"], "GENIC_UNRESOLVED")
        self.assertEqual(self.annotate(VEP_MATCH_STATUS="NO_GENE_MATCH")["DOWNSTREAM_CONTEXT"], "GENIC_UNRESOLVED")

    def test_inversion_gene_span_is_not_breakpoint_disruption(self):
        row = self.annotate(SVTYPE="INV", SV_GENE_RELATIONSHIP="INVERSION_SPANS_INTACT_GENE")
        self.assertEqual(row["DOWNSTREAM_INV_EFFECT_SCOPE"], "INTACT_GENE_SPAN_CONTEXT")
        self.assertEqual(row["DOWNSTREAM_CONTEXT_SCOPE"], "ANNOTATED_FEATURES_NOT_EXACT_BREAKPOINTS")

    def test_inversion_support_uses_maximum_not_sum(self):
        row = self.annotate(SVTYPE="INV", CALLER_READ_SUPPORT="Sniffles2:2;cuteSV:2;Delly:2")
        self.assertEqual(row["DOWNSTREAM_MAX_READ_SUPPORT"], 2)
        self.assertIn("INV_READ_SUPPORT_BELOW_MINIMUM", row["DOWNSTREAM_EXCLUSION_REASONS"])
        for reads, group in [(3, "MINIMUM_TO_STRONG"), (4, "MINIMUM_TO_STRONG"), (5, "AT_OR_ABOVE_STRONG")]:
            row = self.annotate(SVTYPE="INV", EVENT_MAX_CALLER_READ_SUPPORT=str(reads))
            self.assertEqual(row["DOWNSTREAM_INV_READ_GROUP"], group)
        row = self.annotate(SVTYPE="INV", CALLER_READ_SUPPORT=".")
        self.assertIn("INV_READ_SUPPORT_UNKNOWN", row["DOWNSTREAM_REVIEW_REASONS"])

    def test_structured_read_support_and_caller_count_conflicts(self):
        count, source = filtering.max_read_support({"CALLER_EVIDENCE_JSON": json.dumps([{"CALLER_SUPPORT": "3"}, {"CALLER_SUPPORT": "5"}])})
        self.assertEqual((count, source), (5, "CALLER_EVIDENCE_JSON"))
        self.assertIn("CALLER_COUNT_CONFLICT", self.annotate(CALLER_COUNT="2")["DOWNSTREAM_REVIEW_REASONS"])

    def test_original_values_ranks_and_bytes_are_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "input.tsv"
            rows = [{**base_row(), "SV_ID": f"event{i}", "FINAL_EVENT_RANK_GLOBAL": str(i), "CALLER_COUNT": str(count), "CALLERS": ".", "CUSTOM_FIELD": "NA"} for i, count in [(1, 1), (2, 3), (3, 2), (4, 3)]]
            pd.DataFrame(rows).to_csv(source, sep="\t", index=False)
            original = source.read_bytes()
            self.run_filter(source, root / "filtered")
            self.assertEqual(source.read_bytes(), original)
            kept = pd.read_csv(root / "filtered/gt_1kb/retained.tsv", sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual(kept[list(rows[0])].to_dict("records"), rows)
            priority = pd.read_csv(root / "filtered/gt_1kb/caller_priority.tsv", sep="\t", dtype=str)
            self.assertEqual(list(priority["FINAL_EVENT_RANK_GLOBAL"]), ["2", "4", "3", "1"])

    def test_unknown_annotation_and_length_are_reviewed_not_discarded(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            rows = [base_row(), {**base_row(), "SV_ID": "unknown_context", "VEP_TRANSCRIPT_REGION_CLASS": "."},
                    {**base_row(), "SV_ID": "bnd", "SVTYPE": "BND"}]
            pd.DataFrame(rows).to_csv(root / "input.tsv", sep="\t", index=False)
            audit = self.run_filter(root / "input.tsv", root / "filtered").set_index("SV_ID")
            self.assertEqual(audit.loc["unknown_context", "DOWNSTREAM_DECISION"], "REVIEW_REQUIRED")
            self.assertEqual(audit.loc["bnd", "DOWNSTREAM_DECISION"], "REVIEW_REQUIRED")
            self.assertEqual(len(pd.read_csv(root / "filtered/length_unresolved/review_required.tsv", sep="\t")), 1)

    def test_event_wide_conflict_is_not_resolved_by_gene_rank(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            pd.DataFrame([base_row(), {**base_row(), "GENE": "OTHER_GENE", "NEEDLR_AF": "0.005"}]).to_csv(root / "input.tsv", sep="\t", index=False)
            audit = self.run_filter(root / "input.tsv", root / "filtered")
            self.assertTrue(audit["DOWNSTREAM_DECISION"].eq("REVIEW_REQUIRED").all())
            self.assertTrue(audit["DOWNSTREAM_REVIEW_REASONS"].str.contains("EVENT_WIDE_EVIDENCE_CONFLICT").all())

    def test_rerun_removes_only_owned_exports_and_keeps_notes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source, out = root / "input.tsv", root / "filtered"
            pd.DataFrame([base_row()]).to_csv(source, sep="\t", index=False)
            self.run_filter(source, out)
            notes = out / "my_review_notes.txt"
            notes.write_text("Keep these notes")
            self.assertTrue((out / "gt_1kb/by_sv_type/DEL/all.tsv").exists())
            pd.DataFrame([{**base_row(), "SVTYPE": "DUP", "SVLEN": "2000"}]).to_csv(source, sep="\t", index=False)
            self.run_filter(source, out)
            self.assertFalse((out / "gt_1kb/by_sv_type/DEL/all.tsv").exists())
            self.assertEqual(notes.read_text(), "Keep these notes")

    def test_actual_compact_table_and_separate_small_context_files(self):
        from build_candidate_tables import build_sv_table
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            rows = [{**base_row(), "GENES": "TEST_GENE", "SVLEN": "-999", "SV_ID": "small_exon"},
                    {**base_row(), "GENES": "OTHER_GENE", "GENE": "OTHER_GENE", "SVLEN": "-999", "SV_ID": "small_intron", "VEP_TRANSCRIPT_REGION_CLASS": "INTRONIC"}]
            compact = build_sv_table(pd.DataFrame(rows), 10000)
            compact.to_csv(root / "input.tsv", sep="\t", index=False)
            self.run_filter(root / "input.tsv", root / "filtered")
            self.assertTrue((root / "filtered/lt_1kb/by_sv_type/DEL/exonic_or_splice.tsv").exists())
            self.assertTrue((root / "filtered/lt_1kb/by_sv_type/DEL/intronic.tsv").exists())

    def test_empty_and_populated_runs_produce_valid_separate_figures(self):
        for empty in [True, False]:
            with self.subTest(empty=empty), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                data = pd.DataFrame([base_row()])
                if empty:
                    data = data.iloc[:0]
                else:
                    data = pd.concat([data, pd.DataFrame([{**base_row(), "SV_ID": "small", "SVLEN": "-999", "NEEDLR_AF": "."}])], ignore_index=True)
                data.to_csv(root / "input.tsv", sep="\t", index=False)
                self.run_filter(root / "input.tsv", root / "filtered")
                save = plotting.save_figure
                with patch.object(sys, "argv", ["plot", "--filter-dir", str(root / "filtered"), "--out-dir", str(root / "plots")]), patch.object(plotting, "save_figure", side_effect=lambda fig, path, **kw: save(fig, path, dpi=60)), contextlib.redirect_stdout(io.StringIO()):
                    plotting.main()
                self.assertEqual(len(list((root / "plots").rglob("*.png"))), 4)
                for png in (root / "plots").rglob("*.png"):
                    with Image.open(png) as image:
                        image.verify()
                    ET.parse(png.with_suffix(".svg"))
                    self.assertTrue(png.with_suffix(".pdf").read_bytes().rstrip().endswith(b"%%EOF"))


if __name__ == "__main__":
    unittest.main()
