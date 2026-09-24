# Thesis plotting layer

This directory contains downstream visualization and conservative intersection scripts for the SV-PIPELINE outputs. The plotting code is intentionally separated from variant calling and annotation so that figures can be regenerated without rerunning the biological pipeline.

## Design principles

- `SV_ID` is the unit for variant-burden plots. The integrated LRS/SRS table contains one row per `(SV_ID, overlapping gene)`, so burden plots deduplicate master SVs.
- `(SV_ID, gene)` is retained for candidate-prioritization and evidence-integration plots.
- needLR remains a supplementary population-frequency branch rather than the master SV universe.
- Straglr, TLDR, methylation and phasing are complementary biological layers and are not forced into the Jasmine caller-concordance calculation.
- Straglr/TLDR can be attached to the master table only through explicit coordinate-aware matching with `intersect_orthogonal_sv_evidence.py`.
- Every figure is exported as PDF, SVG and 600-dpi PNG.
- Plot scripts also write the summarized TSV used to make the figure when appropriate.

## Scripts

- `plot_lrs_qc.py`: mosdepth coverage QC.
- `plot_caller_concordance.py`: caller-aware UpSet-style plot for Jasmine or SURVIVOR summaries.
- `plot_sv_landscape.py`: SV type, size, chromosome distribution and caller-support landscape.
- `plot_needlr_population.py`: overall and ancestry-specific needLR population-frequency plots.
- `plot_candidate_genes.py`: gene-prioritization plot from `*_ranked_candidates.tsv`.
- `plot_gene_hpo_heatmap.py`: human gene-HPO association heatmap from the offline Monarch branch.
- `plot_candidate_evidence_matrix.py`: integrated SV/gene evidence matrix from the master TSV.
- `plot_straglr.py`: tandem-repeat locus size/copy-number/support overview.
- `plot_mei.py`: TLDR mobile-element insertion summary.
- `plot_methylation.py`: candidate-region modkit methylation track; bedMethyl from `modkit pileup` is preferred for final thesis figures.
- `plot_phasing_qc.py`: WhatsHap/LongPhase phased-genotype QC.
- `plot_lrs_vs_srs.py`: Truvari-based LRS-vs-SRS concordance visualization.
- `intersect_orthogonal_sv_evidence.py`: conservative Straglr/TLDR coordinate matching to the master integrated table.
- `run_thesis_plots.py`: launch all applicable LRS or SRS plots for one sample using the current pipeline directory structure.
- `plot_utils.py`: shared styling, parsing and figure-export helpers.

## Current status

These scripts are written against the output contracts currently defined in the repository. Real output files are not yet available, so column aliases are handled conservatively and clear errors are raised when a required field cannot be inferred. Once real pipeline outputs are available, adjust column mappings where required rather than changing the biological logic.

The plotting layer is intentionally not added to `rule all` yet. It should remain downstream until the first real outputs have been inspected and the final column contracts confirmed.

## Environment

```bash
conda env create -f envs/plots.yaml
conda activate svplots
```

## Run all available plots for one LRS sample

```bash
python plots/run_thesis_plots.py \
  --root /path/to/outs \
  --sample patient01
```

For an SRS sample, pass the platform so `SUPP_VEC` is decoded as Manta,DELLY
and LRS-only plots are skipped:

```bash
python plots/run_thesis_plots.py \
  --root /path/to/srs_results \
  --sample patient01 \
  --platform srs
```

Optional candidate-region methylation:

```bash
python plots/run_thesis_plots.py \
  --root /path/to/outs \
  --sample patient01 \
  --methylation-region chr3:193600000-193670000
```

Dry-run discovery without generating figures:

```bash
python plots/run_thesis_plots.py \
  --root /path/to/outs \
  --sample patient01 \
  --dry-run
```

## Individual examples

```bash
python plots/plot_caller_concordance.py \
  --input outs/patient01/sv/merged/patient01_caller_support_summary.tsv \
  --out-prefix outs/patient01/plots/patient01_caller_concordance

python plots/plot_sv_landscape.py \
  --input outs/patient01/gene_discovery/patient01_integrated_SV_gene_analysis.tsv \
  --out-prefix outs/patient01/plots/patient01_sv_landscape
```

## Orthogonal evidence intersection

```bash
python plots/intersect_orthogonal_sv_evidence.py \
  --integrated outs/patient01/gene_discovery/patient01_integrated_SV_gene_analysis.tsv \
  --straglr outs/patient01/sv/straglr/patient01_straglr.annotated.tsv \
  --tldr outs/patient01/mei/tldr/patient01.tldr.table.txt \
  --output outs/patient01/gene_discovery/patient01_integrated_with_orthogonal_evidence.tsv
```

For thesis figures, use PDF or SVG as the primary figure and PNG only when raster output is required.


## Snakemake integration

The core LRS workflow and downstream interpretation are intentionally separated.

Run the biological pipeline with:

```bash
cd snakemake_pipelines/lrs
snakemake -s Snakefile_LRS_update --use-conda --cores 32
```

Then build the orthogonal-evidence interpretation tables with:

```bash
snakemake -s Snakefile_LRS_postprocess --use-conda --cores 8
```

The post-processing `rule all` generates only:

```text
<sample>_integrated_SV_gene_with_orthogonal_evidence.tsv
```

Thesis plots are a separate explicit target and are never part of either core `rule all`:

```bash
snakemake -s Snakefile_LRS_postprocess --use-conda --cores 8 all_thesis_plots
```

For one sample only, target its marker directly:

```bash
snakemake -s Snakefile_LRS_postprocess --use-conda --cores 8 \
  "<output_root>/<sample>/plots/.thesis_plots.done"
```

If a candidate methylation region has been selected, set `thesis_methylation_region`
in `config_lrs.yaml`; otherwise methylation plotting is skipped while the other plots run.


## Organized output structure

The sample-level runner now writes figures and their source TSVs into thematic folders:

```text
<sample>/plots/
├── 01_qc/
├── 02_caller_concordance/
├── 03_sv_landscape/
├── 04_population_frequency/
├── 05_candidate_prioritization/
├── 06_phenotype/
├── 07_orthogonal/
│   ├── straglr/
│   └── tldr/
├── 08_phasing/
├── 09_methylation/
└── 10_integrated_evidence/
```

The post-processing workflow also creates three interpretation tables in
`<sample>/gene_discovery/`:

```text
<sample>_integrated_SV_gene_with_orthogonal_evidence.tsv
<sample>_independent_orthogonal_findings.tsv
<sample>_gene_multimodal_evidence_summary.tsv
```

The extended master table keeps the Jasmine-defined SV universe and adds
coordinate-aware Straglr/TLDR evidence plus LongPhase SV phasing. The
independent table retains PASS TLDR and Straglr findings without a Jasmine
counterpart. The gene summary collapses evidence without double-counting the
same master SV across multiple SV-gene rows.

WhatsHap is retained as a small-variant phasing/QC layer rather than treated as
an SV caller. Methylation remains a candidate-region annotation because a CpG
overlap is not equivalent to SV confirmation.
