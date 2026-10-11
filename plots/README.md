# Thesis plotting layer

This directory contains downstream visualization and conservative intersection scripts for the SV-PIPELINE outputs. The plotting code is intentionally separated from variant calling and annotation so that figures can be regenerated without rerunning the biological pipeline.

## Design principles

- `SV_ID` is the unit for variant-burden plots. The integrated LRS/SRS table contains one row per `(SV_ID, overlapping gene)`, so burden plots deduplicate master SVs.
- `(SV_ID, gene)` is retained for candidate-prioritization and evidence-integration plots.
- needLR remains a supplementary population-frequency branch rather than the master SV universe.
- Straglr, TLDR, methylation and phasing are complementary biological layers and are not forced into the Jasmine caller-concordance calculation.
- Straglr/TLDR can be attached to the master table only through explicit coordinate-aware matching with `intersect_complementary_sv_evidence.py`.
- Figures are exported as PDF, SVG and PNG: normally 600 dpi for overviews,
  300 dpi for loci and 200 dpi for large candidate cards.
- Plot scripts also write the summarized TSV used to make the figure when appropriate.

## Scripts

- `plot_lrs_qc.py`: mosdepth coverage QC.
- `plot_caller_concordance.py`: caller-aware UpSet-style plot for Jasmine or SURVIVOR summaries.
- `plot_sv_landscape.py`: SV type, size, chromosome distribution and caller-support landscape.
- `plot_genome_wide_sv_map.py`: chromosome position strips, BND/TRA arcs,
  rainfall spacing, optional measured depth and repeat comparisons.
- `plot_needlr_population.py`: overall and ancestry-specific needLR population-frequency plots.
- `plot_candidate_genes.py`: final gene-prioritization figures written separately for panel and non-panel genes; inheritance class and best SV mechanism are shown explicitly when available.
- `plot_gene_hpo_heatmap.py`: human gene-HPO association heatmap, plus a separate teal/purple custom-versus-Exomiser phenotype-rank comparison. The comparison uses the same HPO query and comparable gene set; it does not replace the integrated ranking. Its source TSV retains unavailable genes and marks which genes are displayed.
- `plot_candidate_evidence_matrix.py`: integrated SV/gene evidence matrix from the compact final SV-gene table, selected by final panel-specific rank and displaying AD/AR/mechanism evidence.
- `plot_panel_nonpanel_discovery.py`: explicit panel versus non-panel SV-gene categories, separating HPO and curated disease evidence.
- `plot_sv_gene_associations.py`: quantitative relationship between technical SV support, phenotype relevance, gene-disease evidence and each explicit SV-gene pair.
- `plot_sv_gene_network.py`: bipartite map showing which master SV overlaps which gene, including multi-gene SVs and genes affected by multiple prioritized SVs.
- `plot_candidate_locus.py`: one detailed locus figure per prioritized SV-gene pair, including nearby genes, principal evidence fields and optional indexed modkit methylation.
- `run_samplot_candidates.py`: read-level ONT breakpoint/alignment evidence for the highest-priority master SVs; large events are breakpoint-zoomed and BNDs are shown as breakpoint-context views.
- `plot_mitocarta_sv_genes.py`: MitoCarta3.0 view of nuclear-encoded mitochondrial genes intersected by SVs, including pathway context separated by direct/proximal/interval SV-gene relationship.
- `plot_mitochondrial_gene_ranking.py`: separate ranked figures for nuclear-encoded mitochondrial genes and mtDNA-encoded genes.
- `plot_cohort_sv_comparison.py`: deduplicated comparison of completed samples by SV type and evidence-defined subsets.
- `plot_straglr.py`: tandem-repeat locus size/copy-number/support overview.
- `plot_mei.py`: TLDR mobile-element insertion summary.
- `plot_methylation.py`: candidate-region modkit methylation track; bedMethyl from `modkit pileup` is preferred for final thesis figures.
- `plot_phasing_qc.py`: WhatsHap/LongPhase phased-genotype QC.
- `plot_lrs_vs_srs.py`: Truvari-based LRS-vs-SRS concordance visualization.
- `intersect_complementary_sv_evidence.py`: conservative Straglr/TLDR coordinate matching and SV-aware LongPhase evidence. This is the active LRS implementation; the older duplicate `intersect_orthogonal_sv_evidence.py` has been removed.
- `run_thesis_plots.py`: launch all applicable LRS or SRS plots for one sample using the current pipeline directory structure.
- `plot_utils.py`: shared styling, parsing and figure-export helpers.

## Current status

### Added genome and candidate views

Only one new production script was needed; other additions extend the existing
scripts and launcher. Both LRS and SRS use the shared plotting code. Source TSVs
accompany the figures.

| View | Folder within a sample's `plots/` | How to read it |
|---|---|---|
| Position map and rainfall | `03_sv_landscape/` | Same Mb scale for every chromosome; SV type colours; BND/TRA arcs use the reported partner. Rainfall uses one START anchor per master SV and restarts at each chromosome. |
| Counts and evidence UpSet | `00_start_here/analysis_counts/` | Raw caller records and unique merged SVs have separate panels. Evidence groups overlap; they are not sequential filters. |
| Size versus support | `03_sv_landscape/` | Black rings highlight shortlisted SVs. Caller counts and read support are separate; read counts are not summed across callers. |
| Population AF with candidates | `04_population_frequency/` | Measured AF <=1% and unknown AF stay in one review group, with different colours. Histograms exclude zero/missing values and state their counts; black marks identify candidates. |
| HP1/HP2 gene view | `12_candidate_loci/top_N/gene_phase/` | Only compatible phased SNV/SV pairs on the same chromosome and phase block are placed on haplotypes. Other pairs remain unresolved. |
| One-page candidate cards | `12_candidate_loci/top_N/candidate_cards/` | A2 portrait PDFs combine available locus, depth, methylation, Samplot, phase and optional assembly evidence. Originals remain available for detail. |
| Repeat context | `03_sv_landscape/*repeat_context*` | Observed breakpoint overlap versus matched random positions; needs repeat BED files. |
| Local assembly alignment | `12_candidate_loci/top_N/*assembly*` | Forward/reverse PAF blocks use different colours; needs an existing local assembly alignment. |

The map includes canonical chromosomes present in the actual reference index,
including X/Y/MT. Other contigs, invalid positions and conflicting propagated
coordinates remain in the source TSV with their status. This is separate from
the nuclear filtering views. Genotype/VAF, when present, are retained in the
size/support TSV; none is inferred from population AF or caller counts.

**Interpretation:** a cluster can reflect mapping difficulty or biology. Repeat
overlap does not establish NAHR or another mechanism. Missing AF does not
establish a novel SV. A shortlist is a review workload, not proven diagnoses.
Patient HPO counts require explicit patient-specific scope; HON seed terms do
not qualify. Unknown evidence is reported separately, not assigned support.

### Measured depth and optional reference tracks

Snakemake passes the actual FASTA index, panel BED and DELLY exclusion BED.
LRS reuses its existing binned-depth rule. SRS saves 10 kb mean-depth bins in
its existing mosdepth QC step, using `--by` and skipping the large per-base
file. Its read-selection defaults and summary filenames are preserved. Older
SRS results need that QC step rerun once for the new `.regions.bed.gz` output.

The genome ribbon pools measured bins into 1 Mb bins for visibility. Autosomes
share a length-weighted autosomal mean baseline; X/Y/MT use their own chromosome
mean. A ratio of 1 on X/Y therefore does not establish expected chromosome copy
number. Bins touching supplied exclusions/gaps/centromeres are masked. The TSV
records the baseline and covered bases. This is relative depth, not absolute
copy number or GC-corrected CNVpytor. For individual DEL/DUP events, inspect
flank-normalised local depth (LRS), local measured mean depth (SRS) and Samplot;
the 1 Mb ribbon cannot resolve small events. Missing bins are never invented.

Both configs accept matching GRCh38 tracks:

```yaml
genome_plot_centromere_bed: /path/to/hg38_centromeres.bed
genome_plot_gap_bed: /path/to/hg38_gaps.bed
genome_plot_repeat_beds:
  Tandem_repeats: /path/to/SV-PIPELINE/reference/human_GRCh38_no_alt_analysis_set.trf.bed
  Segmental_duplications: /path/to/hg38_segdup.bed
  LINE: /path/to/hg38_LINE.bed
  SINE: /path/to/hg38_SINE.bed
genome_plot_random_draws: 200
```

Defaults are null/empty: unavailable tracks are labelled, and repeat analysis
is skipped. Download resources manually with the
[UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables): choose Human,
**GRCh38/hg38**, whole genome, and BED output. Use `gap` for assembly gaps;
`cytoBand`, filtered to `gieStain = acen`, for cytogenetic centromere context;
`genomicSuperDups` for segmental duplications; and separate `rmsk` exports
filtered by `repClass = LINE` or `SINE`. Cytobands are not precise molecular
centromere boundaries. Export BED, not the full UCSC database table. Check the
assembly against your aligned reference: bounds checks cannot prove assembly
identity. BED is 0-based/end-excluded; master START/POS2 are VCF 1-based.

Repeat comparisons use 100 bp windows around breakpoints and count each event
once per class if either endpoint overlaps. Random positions retain chromosome,
SV type and interval span; intrachromosomal BND separation is retained, and
interchromosomal BND endpoints remain on their respective chromosomes. Both
observed and random endpoints avoid the same supplied exclusion mask;
unresolved BND partners are excluded. Black marks show the mean and central
95% of the random draws (seed 27), not a confidence interval for a mechanism
or an adjusted significance test. Classes can overlap. No read-callability or
mappability mask is inferred, so this comparison is exploratory. Several repeat
classes for many SVs may take minutes. TLDR/MELT findings remain separate from
repeat-overlap percentages and retain their existing candidate evidence states.

### Phasing and optional local assembly

LRS gene phase reuses `*_snv_sv_candidates.tsv` and joint LongPhase SV genotypes.
HP1/HP2 are local block labels, not maternal/paternal labels. Pairing now checks
chromosome identity and refuses to borrow phase from another ALT at the same
position. Multiallelic records stay unresolved. SRS currently phases small
variants only; those calls cannot place an SV on HP1/HP2. Haplotag-coloured
Samplot is not enabled: compatible SV-aware HP tags in the BAM would be needed
before that display could support a cis/trans interpretation.

For key candidates with an **existing local assembly**, align its contigs to
the same reference. If minimap2 is installed:

```bash
minimap2 -cx asm5 /path/to/matching_reference.fa candidate_assembly.fa > candidate.paf
```

See [official minimap2 instructions](https://github.com/lh3/minimap2) for
installation and alignment choices. Assembly generation is a separate analysis
decision; no assembly/download/install bootstrap script has been added. Create
a TSV with one `SV_ID` and absolute `PAF_PATH` per candidate and set
`candidate_assembly_manifest: /path/to/assemblies.tsv` in the relevant config.
Plots join PAF block endpoints on separate contig/reference-chromosome axes;
these are alignment-block dotplots, not base-by-base identity matrices or a
validated reconstruction.

For completed results, the existing launcher accepts the same inputs:

```bash
python plots/run_thesis_plots.py --root /path/to/results --sample patient01 \
  --platform srs --gene-bed /path/to/genes.bed \
  --reference-fai /path/to/reference.fa.fai --panel-bed /path/to/panel.bed \
  --exclude-bed /path/to/delly_hg38_exclude.bed
```

Use `--platform lrs` for LRS. Optional flags: `--centromere-bed`, `--gap-bed`,
`--depth-bed`, repeated `--repeat-bed CLASS=PATH`, and `--assembly-manifest`.
Samplot images are attached through their manifest's exact SV IDs. Snakemake
waits for the current Samplot manifest before making cards. If Samplot is
created separately later, rerun the launcher to refresh cards.

Technical references: [BED coordinates](https://genome.ucsc.edu/FAQ/FAQformat),
[mosdepth bins](https://github.com/brentp/mosdepth), and
[phase sets](https://whatshap.readthedocs.io/en/latest/guide.html).

The plotting layer is written against the current integrated LRS/SRS output contracts and is now being used with completed LRS sample outputs. Column aliases remain conservative so schema changes fail visibly rather than silently changing biological interpretation.

The plotting layer remains intentionally downstream and is not part of the biological pipeline `rule all`. Figures can therefore be regenerated or extended without rerunning SV calling or annotation.

## Environment

For a small initial review of completed outputs, start with
[Interpreting LRS results](../docs/INTERPRETING_LRS_RESULTS.md). The launcher now
supports `--review-only`, creates `00_start_here/START_HERE.txt` and an editable
review TSV, separates panel/non-panel evidence matrices, and adds evidence
availability and mechanism summaries. Optional `--gene-gtf` provides actual
exon tracks; existing depth bins are used where available. Existing review
notes are preserved on rerun.

After optional further filtering, `plot_downstream_sv_filtering.py` creates
four separate figures explaining filter outcomes, large/small annotation
context and unique-SV caller/AF groups. See
[Downstream SV filtering](../docs/DOWNSTREAM_SV_FILTERING.md) for the commands,
counting units and scientific limits.

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

## Complementary evidence intersection

```bash
python plots/intersect_complementary_sv_evidence.py \
  --integrated outs/patient01/gene_discovery/patient01_integrated_SV_gene_analysis.tsv \
  --straglr outs/patient01/sv/straglr/patient01_straglr.annotated.tsv \
  --tldr outs/patient01/mei/tldr/patient01.tldr.table.txt \
  --output outs/patient01/gene_discovery/patient01_integrated_with_complementary_evidence.tsv
```

For thesis figures, use PDF or SVG as the primary figure and PNG only when raster output is required.


## Snakemake integration

The core LRS workflow and downstream interpretation are intentionally separated.

Run the biological pipeline with:

```bash
cd snakemake_pipelines/lrs
snakemake -s Snakefile_LRS_update --use-conda --cores 32
```

Then build the complementary-evidence interpretation tables with:

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

