"""Scientific reporting regressions; every event in these fixtures is synthetic."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "plots"))
sys.path.insert(0, str(ROOT / "scripts"))
import candidate_review as review
import build_candidate_review as builder
import plot_candidate_evidence_matrix as matrix
import plot_event_review_summary as summary
import plot_candidate_locus as locus
from build_candidate_tables import build_sv_table


def synthetic_rows():
    base = {"SV_ID": "del1", "GENES": "TEST_A", "CHROM": "chr1",
            "START": "100000", "END": "210000", "SVTYPE": "DEL",
            "SVLEN": "-110000", "SV_EVENT_SPAN_BP": "110000",
            "PANEL_STATUS": "PANEL_GENE", "CALLER_COUNT": "2",
            "CALLERS": "Sniffles2;cuteSV", "EVENT_MAX_CALLER_READ_SUPPORT": "18",
            "GENE_RELEVANCE_TIER": "HIGH", "GENE_RELEVANCE_DISPLAY_SCORE": "12",
            "SV_GENE_EFFECT": "WHOLE_GENE_DELETION",
            "INHERITANCE_MECHANISM_CLASS": "SUPPORTED_DISEASE_MECHANISM",
            "EVENT_RANK_WITHIN_PANEL_STATUS": "1", "NEEDLR_AF": "0",
            "NEEDLR_STATUS": "MATCH", "GNOMAD_SV_EXACT_MATCH": "YES",
            "GNOMAD_SV_AF": "0.08", "STRAGLR_MATCH": "NO",
            "TLDR_MATCH": "NOT_AVAILABLE", "LONGPHASE_GT": "0|1",
            "LONGPHASE_PS": "100", "LONGPHASE_PHASED": "YES",
            "DEPTH_PATTERN": "CONSISTENT_WITH_LOSS", "DEPTH_RATIO": "0.55",
            "DEPTH_EVALUATION": "EVALUATED", "METHYLATION_CONTEXT": "EVALUATED"}
    return [base, {**base, "SV_ID": "inv1", "GENES": "TEST_B", "SVTYPE": "INV",
                   "SV_GENE_EFFECT": "INVERSION_SPANS_INTACT_GENE", "EVENT_RANK_WITHIN_PANEL_STATUS": "2",
                   "INHERITANCE_MECHANISM_CLASS": "INTERVAL_CONTEXT_ONLY"},
            {**base, "SV_ID": "dup1", "GENES": "TEST_C", "SVTYPE": "DUP",
             "PANEL_STATUS": "NONPANEL_GENE", "SV_GENE_EFFECT": "PARTIAL_GENE_DUPLICATION",
             "NEEDLR_AF": ".", "NEEDLR_STATUS": "NO_MATCH", "GNOMAD_SV_AF": ".",
             "DEPTH_PATTERN": "CONSISTENT_WITH_LOSS", "RECESSIVE_PAIR_STATUS": "AR_CIS_NOT_BIALLELIC_BY_PHASE"}]


class CandidateReviewTests(unittest.TestCase):
    def test_evidence_states_do_not_invent_negative_results(self):
        empty = review.event_evidence({"SVTYPE": "DEL"})
        self.assertEqual(empty["Straglr"][0], "UNAVAILABLE")
        self.assertEqual(empty["SV phase"][0], "UNAVAILABLE")
        self.assertEqual(review.event_evidence({"STRAGLR": "NO"})["Straglr"][0], "NO_MATCH")
        self.assertEqual(review.event_evidence({"SVTYPE": "INV"})["Local depth"][0], "NOT_APPLICABLE")
        self.assertEqual(review.event_evidence({})["Local depth"][0], "UNAVAILABLE")
        self.assertEqual(review.event_evidence({"LONGPHASE_MATCH": "YES", "LONGPHASE_GT": "0/1"})["SV phase"][0], "UNRESOLVED")
        self.assertEqual(review.event_evidence({"LONGPHASE_MATCH": "NO", "LONGPHASE_GT": "0|1", "LONGPHASE_PS": "100"})["SV phase"][0], "UNRESOLVED")
        self.assertEqual(review.event_evidence({"CALLER_EVIDENCE_FLAGS": "LOW_GQ"})["Call support"][0], "REVIEW")

    def test_population_uses_max_valid_af_and_preserves_zero(self):
        self.assertEqual(review.population({"NEEDLR_AF": "0"}), ("REPORTED", "max AF\n0"))
        self.assertEqual(review.population({"NEEDLR_AF": "0", "GNOMAD_SV_AF": "0.02;0.1"}), ("REVIEW", "max AF\n0.1"))
        self.assertEqual(review.population({"NEEDLR_AF": "0.7", "NEEDLR_STATUS": "NO_MATCH"})[0], "NO_MATCH")
        self.assertEqual(review.population({"NEEDLR_AF": "0.7", "NEEDLR_STATUS": "NOT_EVALUABLE_GE10MB"})[0], "UNAVAILABLE")
        self.assertEqual(review.population({"GNOMAD_SV_AF": "0.7", "GNOMAD_SV_EXACT_MATCH": "NO"})[0], "UNAVAILABLE")

    def test_depth_direction_and_pair_phase_remain_distinct(self):
        obs = review.event_evidence(synthetic_rows()[2])
        self.assertEqual(obs["Local depth"][0], "REVIEW")
        self.assertEqual(obs["SNV + SV"], ("REVIEW", "Cis pair"))
        trans = review.event_evidence({"RECESSIVE_PAIR_STATUS": "AR_TRANS_SNV_SV_CANDIDATE"})
        self.assertEqual(trans["SNV + SV"][0], "SUPPORTING")

    def test_existing_rank_separate_panel_and_exact_event_keys(self):
        rows = synthetic_rows()
        rows += [{**rows[0], "SV_ID": "later", "EVENT_RANK_WITHIN_PANEL_STATUS": "9", "GENE_RELEVANCE_DISPLAY_SCORE": "999"},
                 {**rows[0], "SV_ID": "unknown", "PANEL_STATUS": "."}, rows[0]]
        result = review.select_events(pd.DataFrame(rows), 1)
        self.assertEqual(set(result["SV_ID"]), {"del1", "dup1", "unknown"})
        self.assertEqual(result.set_index("SV_ID").loc["unknown", "_panel_group"], "UNRESOLVED_PANEL")
        self.assertEqual(result.set_index("SV_ID").loc["del1", "PLOT_ORDER_BASIS"], "EVENT_RANK_WITHIN_PANEL_STATUS")

    def test_enrichment_rejects_gene_only_and_ambiguous_matches(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "detail.tsv"
            detail = pd.DataFrame([
                {"SV_ID": "a", "GENES": "A", "READ_SUPPORT": "10", "FINAL_GENE_RELEVANCE_TIER": "UPSTREAM"},
                {"SV_ID": "b", "GENES": "A", "READ_SUPPORT": "20"},
                {"SV_ID": "c", "GENES": "C", "READ_SUPPORT": "1"},
                {"SV_ID": "c", "GENES": "C", "READ_SUPPORT": "2"},
            ])
            detail.to_csv(path, sep="\t", index=False)
            source = pd.DataFrame([{"SV_ID": "a", "GENE": "A", "FINAL_GENE_RELEVANCE_TIER": "FINAL"},
                                   {"SV_ID": "b", "GENE": "B"}, {"SV_ID": "c", "GENE": "C"}])
            enriched = review.enrich_events(source, path).set_index("SV_ID")
            self.assertEqual(enriched.loc["a", "READ_SUPPORT"], "10")
            self.assertEqual(enriched.loc["a", "FINAL_GENE_RELEVANCE_TIER"], "FINAL")
            self.assertEqual(enriched.loc["b", "EVIDENCE_JOIN_STATUS"], "NO_EXACT_SV_GENE_ROW")
            self.assertEqual(enriched.loc["c", "EVIDENCE_JOIN_STATUS"], "AMBIGUOUS_SV_GENE_ROWS")
            self.assertTrue(pd.isna(enriched.loc["c", "READ_SUPPORT"]))

    def test_event_wide_counts_do_not_multiply_large_event_genes(self):
        base = synthetic_rows()[0]
        events = review.normalise_events(pd.DataFrame([base, {**base, "GENES": "OTHER_GENE"}]))
        counts = summary.availability_counts(review.evidence_table(events))
        self.assertEqual(set(counts[counts["DOMAIN"].eq("Call support")]["TOTAL"]), {1})
        self.assertEqual(set(counts[counts["DOMAIN"].eq("Gene effect")]["TOTAL"]), {2})
        conflict = events.copy()
        conflict.loc[1, "STRAGLR_MATCH"] = "NOT_AVAILABLE"
        counts = summary.availability_counts(review.evidence_table(conflict))
        row = counts[counts["DOMAIN"].eq("Straglr") & counts["STATE"].eq("UNRESOLVED")].iloc[0]
        self.assertEqual(row["COUNT"], 1)

    def test_geometric_mechanisms_do_not_claim_coding_loss(self):
        self.assertEqual(summary.mechanism_category({"SVTYPE": "DEL", "SV_GENE_EFFECT": "PARTIAL_GENE_DELETION"}), "Partial gene loss")
        self.assertEqual(summary.mechanism_category({"SVTYPE": "INV", "SV_GENE_EFFECT": "INVERSION_SPANS_INTACT_GENE"}), "Gene spanned by inversion")

    def test_review_rerun_preserves_notes(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            events = review.select_events(pd.DataFrame(synthetic_rows()), 5)
            builder.write_review(events, "SYNTHETIC", out, out / "input.tsv", None)
            editable = out / "SYNTHETIC_candidate_review.tsv"
            original = editable.read_text() + "MY NOTES MUST REMAIN\n"
            editable.write_text(original)
            builder.write_review(events.head(1), "SYNTHETIC", out, out / "input.tsv", None)
            self.assertEqual(editable.read_text(), original)
            self.assertEqual(len(pd.read_csv(out / "SYNTHETIC_candidate_review.generated.tsv", sep="\t")), 1)

    def test_real_exons_and_transcript_selection(self):
        with tempfile.TemporaryDirectory() as folder:
            gtf = Path(folder) / "synthetic.gtf"
            gtf.write_text('chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_name "A"; transcript_id "short";\n'
                           'chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_name "A"; transcript_id "long";\n')
            exons = locus.load_selected_exons(gtf, {"A"})
            self.assertEqual(set(exons["start"]), {100})
            chosen, transcript, basis = locus.representative_exons(exons, "A", "chr1", "short")
            self.assertEqual(transcript, "short")
            self.assertTrue(basis.startswith("VEP_CANONICAL"))
            self.assertEqual(locus.representative_exons(exons, "A", "chr1")[1], "long")
            self.assertTrue(locus.load_selected_exons(None, {"A"}).empty)

    def test_methylation_track_does_not_mix_5hmc_with_5mc(self):
        with tempfile.TemporaryDirectory() as folder:
            bed = Path(folder) / "synthetic.bed.gz"
            bed.touch()
            Path(str(bed) + ".tbi").touch()
            lines = ["chr1\t100\t101\tm\t0\t+\t100\t101\t0\t10\t0.25",
                     "chr1\t102\t103\th\t0\t+\t102\t103\t0\t10\t0.50"]
            tabix = SimpleNamespace(fetch=lambda *args: iter(lines), close=lambda: None)
            with patch.dict(sys.modules, {"pysam": SimpleNamespace(TabixFile=lambda *args: tabix)}):
                values = locus.methylation_records(bed, "chr1", 90, 110, "fraction")
            self.assertEqual(len(values), 1)
            self.assertEqual(values[0][1], 25)

    def test_locus_uses_exact_shortlist_and_event_depth_keys(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            rows = synthetic_rows()
            pd.DataFrame(rows).to_csv(out / "input.tsv", sep="\t", index=False)
            pd.DataFrame([{"SV_ID": "del1", "GENE": "TEST_A"}]).to_csv(out / "selected.tsv", sep="\t", index=False)
            (out / "genes.bed").write_text("chr1\t100000\t210000\tTEST_A\nchr1\t200000\t300000\tTEST_B\n")
            pd.DataFrame([{"SV_ID": sv, "CHROM": "chr1", "PLOT_START": "120000", "PLOT_END": "130000", "NORMALIZED_DEPTH": "0.5"} for sv in ["del1", "different_event"]]).to_csv(out / "depth.tsv", sep="\t", index=False)
            save = locus.save_figure
            with patch.object(sys, "argv", ["locus", "--input", str(out / "input.tsv"), "--selection-table", str(out / "selected.tsv"), "--gene-bed", str(out / "genes.bed"), "--depth-bins", str(out / "depth.tsv"), "--out-dir", str(out / "figures")]), patch.object(locus, "save_figure", side_effect=lambda fig, path, **kw: save(fig, path, dpi=60)):
                locus.main()
            manifest = pd.read_csv(out / "figures" / "candidate_locus_manifest.tsv", sep="\t")
            self.assertEqual(list(manifest["SV_ID"]), ["del1"])
            self.assertEqual(manifest.iloc[0]["depth_bins_plotted"], 1)
            self.assertEqual(manifest.iloc[0]["exon_annotation_scope"], "GENE_INTERVAL_ONLY")

    def test_render_split_groups_and_empty_tables(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            for name, data in [("populated", pd.DataFrame(synthetic_rows())), ("empty", pd.DataFrame(columns=["SV_ID", "GENES", "SVTYPE"]))]:
                source = out / f"{name}.tsv"
                data.to_csv(source, sep="\t", index=False)
                prefix = out / name
                save = matrix.save_figure
                with patch.object(sys, "argv", ["matrix", "--input", str(source), "--out-prefix", str(prefix)]), patch.object(matrix, "save_figure", side_effect=lambda fig, path, **kw: save(fig, path, dpi=60)):
                    matrix.main()
                self.assertEqual(len(pd.read_csv(str(prefix) + "_matrix.tsv", sep="\t")), len(data))
                if not data.empty:
                    for group in ["panel", "nonpanel"]:
                        png = Path(str(prefix) + f"_{group}.png")
                        self.assertGreater(png.stat().st_size, 0)
                        with Image.open(png) as image:
                            image.verify()
                        self.assertTrue(Path(str(prefix) + f"_{group}.pdf").read_bytes().rstrip().endswith(b"%%EOF"))
                        ET.parse(Path(str(prefix) + f"_{group}.svg"))
                save = summary.save_figure
                with patch.object(sys, "argv", ["summary", "--input", str(source), "--out-dir", str(out / name), "--sample", "SYNTHETIC"]), patch.object(summary, "save_figure", side_effect=lambda fig, path, **kw: save(fig, path, dpi=60)):
                    summary.main()
                self.assertTrue((out / name / "mechanisms" / "SYNTHETIC_mechanism_summary.tsv").exists())
                with Image.open(out / name / "evidence_availability" / "SYNTHETIC_evidence_availability.png") as image:
                    image.verify()

    def test_actual_compact_schema_is_supported(self):
        compact = build_sv_table(pd.DataFrame(synthetic_rows()), 10000)
        self.assertIn("GENE", compact)
        selected = review.select_events(compact, 5)
        self.assertEqual(set(selected["SV_ID"]), {"del1", "inv1", "dup1"})
        self.assertEqual(review.event_evidence(selected.iloc[0])["Straglr"][0], "NO_MATCH")

    def test_review_only_launcher_does_not_launch_analysis_or_legacy_plots(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            target = out / "SYNTHETIC" / "gene_discovery" / "final"
            target.mkdir(parents=True)
            pd.DataFrame(synthetic_rows()).to_csv(target / "SYNTHETIC_sv_gene_candidates.ranked.tsv", sep="\t", index=False)
            result = subprocess.run([sys.executable, str(ROOT / "plots" / "run_thesis_plots.py"), "--root", str(out), "--sample", "SYNTHETIC", "--review-only", "--dry-run"], text=True, capture_output=True, env={**os.environ, "MPLBACKEND": "Agg"})
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("build_candidate_review.py", result.stdout)
            self.assertIn("plot_event_review_summary.py", result.stdout)
            self.assertNotIn("plot_gene_sv_spectrum.py", result.stdout)
            self.assertNotIn("snakemake", result.stdout)


if __name__ == "__main__":
    unittest.main()
