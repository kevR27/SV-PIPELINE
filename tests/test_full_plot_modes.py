"""Contracts for full plot regeneration without repeating upstream analysis."""
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVED = ROOT / "snakemake_pipelines" / "lrs_downstream_reanalysis"
LRS = ROOT / "snakemake_pipelines" / "lrs"
SRS = ROOT / "snakemake_pipelines" / "srs_wgs_pipeline"


class FullPlotModeTests(unittest.TestCase):
    def test_all_entrypoints_use_full_plot_mode(self):
        paths = [
            ARCHIVED / "Snakefile",
            LRS / "Snakefile_LRS_postprocess",
            SRS / "Snakefile_SRS_WGS",
            ROOT / "plots" / "run_thesis_plots.py",
        ]
        for path in paths:
            with self.subTest(path=path):
                self.assertNotIn("--focused", path.read_text(encoding="utf-8"))

    def test_full_plot_runner_includes_loci_and_mitochondrial_context(self):
        runner = (ROOT / "plots" / "run_thesis_plots.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("plot_candidate_locus.py", runner)
        self.assertIn("plot_mitocarta_sv_genes.py", runner)
        self.assertIn("plot_mitochondrial_gene_ranking.py", runner)
        self.assertIn("plot_breakpoint_gene_tracks.py", runner)
        self.assertIn("SKIPPED_MISSING_INPUT", runner)

    def test_breakpoint_track_count_follows_locus_top_n(self):
        runner = (ROOT / "plots" / "run_thesis_plots.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('"--top-n", str(args.locus_top_n)', runner)
        self.assertNotIn("min(args.locus_top_n, 8)", runner)
        settings = (LRS / "config_lrs.yaml").read_text(encoding="utf-8")
        self.assertIn("candidate_locus_top_n: 50", settings)

    def test_archived_target_requests_all_saved_r_figures(self):
        snakefile = (ARCHIVED / "Snakefile").read_text(encoding="utf-8")
        target = snakefile.split("rule all_saved_thesis_plots:", 1)[1]
        for folder in (
            "plots_r/population_gene_effect",
            "plots_r/candidate_evidence",
            "plots_r/hpo",
            "plots_r/large_sv_depth",
            "plots_r/qc",
            "plots_reanalysis/.thesis_plots.done",
            "plots/15_downstream_filtering/.filter_plots.done",
            "COHORT_PLOT_OUTPUTS",
        ):
            with self.subTest(folder=folder):
                self.assertIn(folder, target)
        self.assertIn("sample=ARCHIVED_QC_SAMPLES", target)

    def test_wrapper_forces_only_plot_producers(self):
        path = ARCHIVED / "run_all_plots.sh"
        subprocess.run(["bash", "-n", str(path)], check=True)
        script = path.read_text(encoding="utf-8")
        self.assertIn("--use-conda", script)
        self.assertIn("--rerun-incomplete", script)
        self.assertIn("all all_saved_thesis_plots", script)
        self.assertIn("--forcerun", script)
        for rule in (
            "archived_thesis_plots_sample",
            "downstream_filtering_plots",
            "r_population_gene_effect",
            "r_candidate_evidence_heatmap",
            "r_hpo_heatmap",
            "r_large_sv_depth",
            "r_sequencing_qc",
        ):
            self.assertIn(rule, script)
        self.assertNotIn(" --forceall", script)
        self.assertNotIn(" --forcerun genomiser", script)


if __name__ == "__main__":
    unittest.main()
