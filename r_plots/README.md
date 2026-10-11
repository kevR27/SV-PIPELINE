# R plots

These plots are an additional publication-style layer. Existing Python figures
under `plots/` remain available.

The R and Python HPO/evidence heatmaps are alternative renderings of the same
analysis, not independent findings. Use one version of each figure in the
thesis. The R layer is retained because it was explicitly requested and has
its own optional target; it is not required for the Exomiser comparison.

The R environment is defined in:

```text
envs/r_plot.yaml
```

Current figures:

- `plot_candidate_evidence_heatmap.R`: ComplexHeatmap view of technical,
  population, panel, mitochondrial, phasing and complementary evidence.
- `plot_hpo_heatmap.R`: ComplexHeatmap gene-HPO matrix with gene relevance
  and panel annotation.
- `plot_large_sv_depth.R`: ggplot2 binned/normalized depth around large
  deletions and duplications.

They are executed through Snakemake by the `final_thesis_analysis` target.
