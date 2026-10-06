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
- `plot_panel_nonpanel_discovery.py`: explicit panel versus non-panel SV-gene categories, separating HPO and curated disease evidence.
- `plot_sv_gene_associations.py`: quantitative relationship between technical SV support, phenotype relevance, gene-disease evidence and each explicit SV-gene pair.
- `plot_sv_gene_network.py`: bipartite map showing which master SV overlaps which gene, including multi-gene SVs and genes affected by multiple prioritized SVs.
- `plot_candidate_locus.py`: one detailed locus figure per prioritized SV-gene pair, including nearby genes, principal evidence fields and optional indexed modkit methylation.
- `run_samplot_candidates.py`: read-level ONT breakpoint/alignment evidence for the highest-priority master SVs; large events are breakpoint-zoomed and BNDs are shown as breakpoint-context views.
- `plot_mitocarta_sv_genes.py`: MitoCarta3.0 view of nuclear-encoded mitochondrial genes intersected by master SVs, including mitochondrial pathways and existing optic-neuropathy phenotype context.
- `plot_cohort_sv_comparison.py`: deduplicated comparison of completed samples by SV type and evidence-defined subsets.
- `plot_straglr.py`: tandem-repeat locus size/copy-number/support overview.
- `plot_mei.py`: TLDR mobile-element insertion summary.
- `plot_methylation.py`: candidate-region modkit methylation track; bedMethyl from `modkit pileup` is preferred for final thesis figures.
- `plot_phasing_qc.py`: WhatsHap/LongPhase phased-genotype QC.
- `plot_lrs_vs_srs.py`: Truvari-based LRS-vs-SRS concordance visualization.
- `intersect_orthogonal_sv_evidence.py`: conservative Straglr/TLDR coordinate matching to the master integrated table.
- `run_thesis_plots.py`: launch all applicable LRS or SRS plots for one sample using the current pipeline directory structure.
- `plot_utils.py`: shared styling, parsing and figure-export helpers.

## Current status

The plotting layer is written against the current integrated LRS/SRS output contracts and is now being used with completed LRS sample outputs. Column aliases remain conservative so schema changes fail visibly rather than silently changing biological interpretation.

The plotting layer remains intentionally downstream and is not part of the biological pipeline `rule all`. Figures can therefore be regenerated or extended without rerunning SV calling or annotation.

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

The post-processing `rule all` retains one authoritative compressed SV-gene table plus small interpretation summaries:

```text
<sample>_integrated_SV_gene_analysis.final.tsv.gz
<sample>_independent_complementary_findings.tsv
<sample>_mitocarta_annotation_summary.tsv
<sample>_gene_multimodal_evidence_summary.tsv
<sample>_mitochondrial_gene_ranking.tsv
```

Large intermediate annotation tables are compact `.delta.tsv.gz` files under `.postprocess_tmp/` and are removed by Snakemake after the final table is built.

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
├── 07_complementary_evidence/
│   ├── straglr/
│   └── tldr/
├── 08_phasing/
├── 09_methylation/
├── 10_integrated_evidence/
├── 11_sv_gene_associations/
├── 12_candidate_loci/
├── 13_breakpoint_evidence/
└── 14_mitochondrial_context/
```

The post-processing workflow keeps one final integrated table and small persistent summaries in `<sample>/gene_discovery/`. Straglr/TLDR/LongPhase, WhatsHap/modkit, MitoCarta and optional exact gnomAD-SV annotations are generated as compact deltas and merged once. This preserves the Jasmine-defined SV set without creating repeated full-table copies.

WhatsHap is retained as a small-variant phasing/QC layer rather than treated as
an SV caller. Methylation remains a candidate-region annotation because a CpG
overlap is not equivalent to SV confirmation.


## Revised organized output layout

The plot runner now writes each analysis layer into a dedicated sample folder:

```text
<sample>/plots/
├── 01_qc/
├── 02_caller_concordance/
├── 03_sv_landscape/
├── 04_population_frequency/
├── 05_candidate_prioritization/
├── 06_phenotype/
├── 07_complementary_evidence/
│   ├── straglr/
│   └── tldr/
├── 08_phasing/
├── 09_methylation/
└── 10_integrated_evidence/
```

The current architecture does not persist serial orthogonal/multimodal/MitoCarta full-table copies. Those annotation layers are compact temporary deltas attached to the same ranked SV-gene rows and merged once into `*_integrated_SV_gene_analysis.final.tsv.gz`.

Modkit plotting now accepts standard bedMethyl content compressed under `.bedmethyl.gz` or generic `.bed.gz` names. With no configured candidate region it produces a genome-wide canonical-chromosome methylation summary; when `thesis_methylation_region` is set, it produces a detailed regional methylation/coverage track instead.


## Detailed SV-gene interpretation

The detailed plotting layer is deliberately keyed to the master `SV_ID`.

`plot_sv_gene_associations.py` keeps one row per `(SV_ID, gene)` and displays
technical support separately from biological prioritization. This makes it
possible to see, for example, whether a highly phenotype-relevant non-panel gene
is carried by a single-caller or multi-caller SV rather than losing the SV-gene
relationship in a gene-only ranking.

`plot_panel_nonpanel_discovery.py` separates panel genes from non-panel genes
with both phenotype and curated disease evidence, phenotype evidence only,
disease evidence only, and other non-panel overlaps. These are discovery
categories, not pathogenicity classes.

`plot_candidate_locus.py` generates separate candidate figures for the top
SV-gene pairs. Each figure shows the master SV coordinates, nearby gene spans,
the overlapping candidate gene, key caller/population/database/phenotype fields,
and indexed modkit methylation points when available. WhatsHap and methylation
remain contextual layers rather than SV confirmations.

When at least two completed post-processing samples are available,
`all_thesis_plots` also creates:

```text
<output_root>/cohort_plots/lrs_sample_comparison.*
```

The comparison deduplicates by `SV_ID` before calculating SV burden and does
not interpret the evidence-defined subsets as sequential filters.


## Samplot read-level breakpoint evidence

`all_thesis_plots` now also generates Samplot images for the highest-priority
master SVs using the original sample BAM. These figures display the long-read
alignment/depth evidence underlying the candidate call. They are visual
technical evidence and are not treated as an independent caller or orthogonal
molecular validation.

For DEL/DUP/INV/INS, the SV interval/type is passed directly to Samplot. Events
at least 1 Mb are breakpoint-zoomed using the configured `samplot_zoom_bp`.
For BND/TRA records, each available breakpoint is plotted as local alignment
context because a single-interval Samplot view cannot by itself represent the
full interchromosomal adjacency.

Outputs:

```text
<sample>/plots/13_breakpoint_evidence/
├── *.png
├── samplot_manifest.tsv
└── .samplot.done
```

## MitoCarta3.0 gene context

MitoCarta annotations are written as a compact temporary delta and merged into:

```text
<sample>_integrated_SV_gene_analysis.final.tsv.gz
```

A persistent diagnostic summary and a separate mitochondrial-gene ranking are also created:

```text
<sample>_mitocarta_annotation_summary.tsv
<sample>_mitochondrial_gene_ranking.tsv
```

The final table contains:

```text
MITOCARTA_STATUS
MITOCARTA_ENCODING
MITOCARTA_DESCRIPTION
MITOCARTA_MAESTRO_SCORE
MITOCARTA_EVIDENCE
MITOCARTA_SUBCOMPARTMENT
MITOCARTA_MITOPATHWAYS
MITOCARTA_TOP_LEVEL_PATHWAYS
MITO_ON_CONTEXT
MITO_ON_ASSOCIATION_CLASS
```

`MITOCARTA_ENCODING=NUCLEAR_MITOCHONDRIAL_GENE` identifies nuclear genes in the full MitoCarta3.0 mitochondrial inventory across mitochondrial pathways and subcompartments. The dedicated mitochondrial ranking additionally intersects chrM SVs with the gene BED so all mtDNA gene classes can be reviewed separately: protein-coding genes, mitochondrial rRNAs and mitochondrial tRNAs. `MITO_ON_CONTEXT=YES` is disease-context co-annotation, not a causality or pathogenicity statement.


## Mechanism-aware SV-gene event ranking

Mechanism-aware event ranking now occurs before the large audit fields are discarded. The ranked table is a temporary compact file under `.postprocess_tmp/`; its relationship/ranking fields are retained in the final integrated table. This preserves one row per `(SV_ID, gene)` association while avoiding another persistent full-table copy. It adds:

```text
SV_GENE_RELATIONSHIP
BREAKPOINT_DISTANCE_TO_GENE_BP
SV_EVENT_SIZE_CLASS
SV_EVENT_SPAN_BP
SV_GENE_COUNT
SV_FOCALITY_CLASS
EVENT_REVIEW_BUCKET
EVENT_INTERPRETATION_SCOPE
BREAKPOINT_DEFINED_EVENT
VERY_LARGE_GE10MB
EVENT_POPULATION_TIER
EVENT_TECHNICAL_TIER
EVENT_GENE_RELEVANCE_SCORE
EVENT_RANK_WITHIN_BUCKET
EVENT_RANK_WITHIN_BUCKET_PANEL_STATUS
```

The event-level gene relevance score is phenotype relevance plus curated
gene-disease evidence. It deliberately excludes the gene-level SV-count
component so an individual event is not rewarded because the same gene happens
to contain several other SV calls.

The principal event buckets are:

```text
SMALL_MEDIUM_CNV_GENE_CANDIDATE
INSERTION_GENE_CANDIDATE
LARGE_CNV_GENE_CANDIDATE
VERY_LARGE_CNV_GENE_CONTEXT
BREAKPOINT_GENE_CANDIDATE
COMPLEX_INTERVAL_CONTEXT
LARGE_COMPLEX_INTERVAL_CONTEXT
BREAKPOINT_EVENT_UNRESOLVED
```

For DEL/DUP/CNV events, transcript overlap is treated as dosage/disruption
context. For INV/BND events, breakpoint overlap is evaluated separately from
genes that merely lie between the breakpoints. An internal gene within a large
inversion is therefore retained as `INTERVAL_CONTEXT_ONLY` and is not treated
as directly disrupted.

The default breakpoint-proximity window is 10 kb and is configured with:

```yaml
sv_gene_breakpoint_tolerance: 10000
```

This threshold is a prioritization aid, not a clinical pathogenicity criterion.

### Large/complex candidate figures

`plot_large_complex_candidates.py` separates gene-directed large
CNV/breakpoint candidates from inversion/BND interval-only context.

`plot_gene_sv_spectrum.py` shows the size spectrum of unique master SVs
contributing to each prioritized gene (<100 kb, 100 kb–1 Mb, 1–10 Mb, ≥10 Mb
and breakends), plus the number of INV/BND events.

Candidate locus and Samplot selection are now stratified across event classes
so small SVs cannot crowd all large and breakpoint-defined events out of
detailed review.
