"""Small synthetic tests for focused, read-free thesis output modules."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def load(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FocusedOutputTests(unittest.TestCase):
    def test_ranking_filter_sensitivity_keeps_unknown_af(self):
        module = load("scripts/check_ranking_stability.py")
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "ranked.tsv"
            pd.DataFrame([
                dict(SV_ID="a", GENE="GENEA", SVTYPE="DEL",
                     FINAL_EVENT_RANK_GLOBAL=1, CALLER_COUNT=2, NEEDLR_AF="0.0001"),
                dict(SV_ID="b", GENE="GENEB", SVTYPE="DUP",
                     FINAL_EVENT_RANK_GLOBAL=2, CALLER_COUNT=1, NEEDLR_AF="."),
                dict(SV_ID="c", GENE="GENEC", SVTYPE="DEL",
                     FINAL_EVENT_RANK_GLOBAL=3, CALLER_COUNT=1, NEEDLR_AF="0.02"),
            ]).to_csv(inp, sep="\t", index=False)
            summary, details, n_events, n_genes = module.analyze(str(inp), 2)
            self.assertEqual((n_events, n_genes), (3, 3))
            base = summary.set_index("scenario")
            self.assertEqual(base.loc["baseline", "shared_top_n"], 2)
            self.assertEqual(base.loc["af_le_0.01_unknown_review", "shared_top_n"], 2)
            self.assertEqual(base.loc["at_least_2_callers", "shared_top_n"], 1)
            self.assertTrue((details["baseline_top_n"] == True).any())

    def test_breakpoint_track_uses_real_bins_and_true_gtf_exons(self):
        module = load("plots/plot_breakpoint_gene_tracks.py")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "final.tsv"
            pd.DataFrame([
                dict(SV_ID="testdel", GENE="GENEA", CHROM="chr3",
                     START=20100, END=50100, SVTYPE="DEL", FINAL_EVENT_RANK_GLOBAL=1)
            ]).to_csv(candidates, sep="\t", index=False)
            sample = module.choose_candidates(candidates, 8)
            self.assertEqual(len(sample), 1)
            depth = root / "depth.bed"
            depth.write_text("".join(
                f"chr3\t{i}\t{i+10000}\t{14 if 20000 <= i < 50000 else 30}\n"
                for i in range(0, 80000, 10000)
            ))
            bed = root / "genes.bed"
            bed.write_text("chr3\t10000\t60000\tGENEA\n")
            gtf = root / "genes.gtf"
            gtf.write_text('chr3\ttest\texon\t22001\t24000\t.\t+\t.\tgene_id "G1"; gene_name "GENEA"; transcript_id "T1";\n')
            depths = module.coverage_tracks(depth,sample)
            genes = module.gene_models(bed,sample)
            exons = module.exon_models(gtf,sample)
            self.assertGreater(len(depths[0]), 3)
            self.assertEqual(genes[0][0][2], "GENEA")
            self.assertEqual(exons[0][0][0], 22000)
            out = root / "fig" / "test"
            module.draw(sample[0], depths[0], genes[0], exons[0], {}, out)
            self.assertTrue(Path(str(out)+".pdf").is_file())
            self.assertTrue(Path(str(out)+".png").is_file())

    def test_cleanup_never_selects_original_sources(self):
        module = load("scripts/cleanup_old_thesis_outputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "RESULTS"
            root.mkdir()
            patient = root / "19545"
            gene = patient / "gene_discovery"
            gene.mkdir(parents=True)
            original = gene / "19545_ranked_candidates.tsv"
            original.write_text("origin\n")
            old = gene / "final"
            old.mkdir()
            (old / "old.tsv").write_text("legacy\n")
            plots = patient / "plots"
            plots.mkdir()
            targets = module.candidates(root,["19545"],"all")
            self.assertIn(old,targets)
            self.assertIn(plots,targets)
            self.assertNotIn(original,targets)
            self.assertNotIn(gene,targets)


if __name__ == "__main__":
    unittest.main()
