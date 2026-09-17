# Thesis plotting layer

This directory contains downstream visualization scripts for the SV-PIPELINE outputs.
The plotting code is intentionally separated from variant calling and annotation so that
figures can be regenerated without rerunning the biological pipeline.

## Design principles

- `SV_ID` is the unit for variant-burden plots. The integrated LRS table contains one
  row per `(SV_ID, overlapping gene)`, so these plots always deduplicate master SVs.
- `(SV_ID, gene)` is retained for candidate-prioritization and evidence-integration plots.
- needLR remains a supplementary population-frequency branch rather than the master
  SV universe.
- Straglr, TLDR, methylation and phasing are treated as complementary biological layers
  and are not forced into the Jasmine caller-concordance calculation.
- Every figure is exported as PDF, SVG and 600-dpi PNG.
- Every script also writes the summarized TSV used to make the figure when appropriate.

## Scripts

- `plot_caller_concordance.py`: three-caller Jasmine UpSet-style plot and caller-support summary.
- `plot_sv_landscape.py`: SV type, size, chromosome distribution and caller-support landscape.
- `plot_needlr_population.py`: overall and ancestry-specific needLR population-frequency plots.
- `plot_candidate_genes.py`: gene-prioritization plot from `*_ranked_candidates.tsv`.
- `plot_candidate_evidence_matrix.py`: integrated SV/gene evidence matrix from the master TSV.
- `plot_straglr.py`: tandem-repeat locus size/copy-number/support overview.
- `plot_mei.py`: TLDR mobile-element insertion summary.
- `plot_methylation.py`: modkit candidate-region methylation plot.
- `plot_lrs_vs_srs.py`: plots an SV-aware LRS-vs-SRS comparison table, ideally produced by Truvari.
- `plot_utils.py`: shared styling and parsing helpers.

## Current status

These scripts are written against the output contracts currently defined in the repository.
Because real output files are not yet available, column aliases are handled conservatively and
clear errors are raised when a required field cannot be inferred. Once real pipeline outputs are
available, update the aliases only where needed rather than changing the biological logic.

## Example

```bash
conda env create -f envs/plots.yaml
conda activate svplots

python plots/plot_caller_concordance.py \
  --input outs/patient01/sv/merged/patient01_caller_support_summary.tsv \
  --out-prefix outs/patient01/plots/patient01_caller_concordance

python plots/plot_sv_landscape.py \
  --input outs/patient01/gene_discovery/patient01_integrated_SV_gene_analysis.tsv \
  --out-prefix outs/patient01/plots/patient01_sv_landscape
```

For thesis figures, use the PDF or SVG output as the primary figure and the PNG only for applications that require raster images.
