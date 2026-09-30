# Final thesis analysis layer

This layer is additive. It does not replace or delete the existing caller,
Jasmine, AnnotSV, needLR, allele-assessment, or Python plotting outputs.

Run it from `snakemake_pipelines/lrs/` with:

```bash
snakemake \
  -s Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8 \
  final_thesis_analysis
```

## Human-readable final tables

The new final folder contains:

```text
<sample>/gene_discovery/final/
├── <sample>_gene_candidates.tsv
├── <sample>_sv_gene_candidates.tsv
├── <sample>_large_sv_depth.tsv
├── <sample>_large_sv_depth_bins.tsv
└── <sample>_snv_sv_candidates.tsv
```

The older detailed files remain unchanged.

### Gene relevance

`GENE_RELEVANCE_SCORE` is deliberately simple:

```text
phenotype relevance + curated gene-disease evidence
```

The number of SVs in a gene is reported, but it does not add points to the new
gene relevance score.

### SV-specific review

The SV-gene table keeps size, call support, population evidence and gene effect
as separate columns. Large events are not rewarded merely because they are
large.

Copy-number-changing SVs are not a separate variant family from SVs. Deletions
and duplications are structural variants that alter copy number. Inversions,
many breakends/translocations and some insertions can be copy-number neutral.

## Probability and statistics

The workflow intentionally does not report a homemade probability that an SV is
pathogenic.

A posterior or conditional pathogenicity probability would require calibrated
priors and feature likelihoods from an adequately sized, independently labeled
pathogenic/benign SV dataset. Caller count, read support, rarity, gene overlap,
phenotype relevance and database evidence are also correlated, which makes a
simple naive-Bayes combination inappropriate.

For the thesis, use:
- transparent evidence fields;
- rankings for candidate triage;
- cohort counts/proportions;
- internal recurrence;
- validated external classifications when available.

Conditional proportions such as "fraction of rare SVs that are multi-caller"
may be reported descriptively, but they are not probabilities that a candidate
causes disease.

If a calibrated phenotype-aware model is desired later, it should be used as a
separate external annotation and validated independently rather than replacing
the transparent tables.

## Large deletion/duplication depth

Mosdepth is run at 10 kb base resolution with MAPQ >=20. Plotting then combines
those base bins into median bins:

- <1 Mb event: 10 kb plot bins
- 1-10 Mb event: 50 kb plot bins
- >=10 Mb event: 100 kb plot bins

This suppresses isolated depth spikes without smoothing the raw BAM itself.

The depth result reports:
- median depth inside the SV;
- median flanking depth;
- inside/flank depth ratio;
- a descriptive loss/gain-consistency flag.

The flag is supporting evidence, not independent molecular confirmation.

## Cohort recurrence

`cohort_analysis/sv_recurrence.tsv` compares highly similar SVs across all
completed study samples.

This is called recurrence, not allele frequency. The study cohort is too small
and selected to be treated as a population-frequency reference.

## Small variant + SV review

The Clair3/LongPhase small-variant branch remains secondary because the thesis
focus is structural variation.

`<sample>_snv_sv_candidates.tsv` reports coding/splice small variants in genes
that also contain a candidate SV. If both variants are phased in the same
LongPhase phase set, the table reports CIS or TRANS. Otherwise it reports
UNRESOLVED.

This is useful for recessive genes but does not replace segregation analysis.

## Nanopore run QC

Dorado QC is optional. Configure a basecalling BAM:

```yaml
dorado_summary_bams:
  "21930": "/path/to/21930.basecalled.bam"
```

The workflow runs:

```bash
dorado summary BASECALL.bam
```

and summarizes yield, read count, read length/N50, Q-score, duration and channel
usage. Do not point this automatically at a later merged BAM unless the required
Dorado tags have been preserved.

## R figures

Publication-style R figures are written in parallel to the existing Python
plots:

```text
<sample>/plots_r/
├── candidate_evidence/
├── hpo/
└── large_sv_depth/
```

The existing `plots/` folder is not removed or overwritten.
