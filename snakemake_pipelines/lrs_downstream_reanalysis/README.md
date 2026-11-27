# LRS downstream reanalysis from saved outputs (without BAM)

This directory runs the **existing LRS downstream interpretation** from
previously generated files. It reuses the rules in
`../lrs/Snakefile_LRS_postprocess`, instead of maintaining a separate,
possibly divergent copy of the scoring algorithms. It **does not rerun**
basecalling, Clair3, SV callers, Jasmine, LongPhase, mosdepth or modkit.
Genomiser is included and mandatory.

**Different from `../genomiser_vcf_only/`:** this workflow rebuilds with retained BAM-derived evidence
 the full saved-data SV gene/event prioritization, multimodal integration,
mitochondrial ranking, cohort recurrence, clinical candidate TSVs, downstream
filtering and selected R/summary plots. It therefore needs more of the
previous run's BAM-derived TSV/VCF/BED output files.

## Selected patients (current study)

The current `config.yaml` explicitly selects `13048`, `19545`,
`11446` and `14996` in this order. The original LRS technical
controls `21930` and `25428` are excluded from reanalysis.
The private patient table may contain other patients; phenotype lookup
uses only the matching sample ID. If the requested sample is actually
`194545`, confirm the ID before changing the mapping, as the clinical
table currently contains `19545`.

## Automatic sample discovery (optional)

Set the common root only once:

```yaml
input_root: "/DATA/casadei7/nanopore/WGS/RESULTS"
sample_ids: auto
```

If `sample_ids: auto` is set explicitly, numeric subfolders (including
control folders such as `21930/` and `25428/`) would be detected automatically.
Do not use auto for the patient-only comparison without excluding controls. Original filenames are derived from the sample directory
name; no BAM paths or individual per-patient VCF paths are needed.

To analyze just one saved patient, change to:

```yaml
sample_ids: ["11446"]
```

For nonnumeric IDs use an explicit list. The workflow inherits the existing
LRS reference/Monarch/MitoCarta/AnnotSV/VEP paths and ranking settings from
`../lrs/config_lrs.yaml`, but it overrides the old sample/BAM mapping.

## Which folders must be restored?

Below `<input_root>/11446/`, restore these **seven directories** and their
actual output files (folders alone are insufficient):

```text
11446/
  snp_clair3/
    11446.vcf.gz
  sv/
    merged/
      11446_merged_SV.vcf.gz
    straglr/
      11446_straglr.annotated.tsv
  gene_discovery/
    11446_integrated_SV_gene_analysis.tsv
    11446_ranked_candidates.tsv
    11446_human_gene_phenotypes.tsv
  mei/
    tldr/
      11446.tldr.table.txt
  phasing_longphase/
    11446.longphase.vcf.gz
    11446.longphase_SV.vcf.gz
  coverage/
    binned/
      11446.10kb.regions.bed.gz
```

The script needs both original genome-wide SNV and SV VCFs. Restore their
indexes (`.tbi` or `.csi`) if available. Additional existing files can
be retained for better multimodal interpretation:
`phasing/<sample>.phased.vcf.gz`,
`methylation/<sample>.cpg.bedmethyl.gz`, caller support summaries,
and any previous QC / saved plots. The workflow regenerates `allele_assessment.tsv` from the saved original
integrated TSV and patient profile. It cannot reconstruct the original integrated
TSV, `ranked_candidates.tsv`,
gene-phenotype tables, or genuine 10 kb read-depth bins from just a
merged SV VCF. Those source files have to be restored from the previous run.

Do not fabricate blank/zero read-depth or phasing files, or treat missing
evidence as negative results.

## Audit BEFORE Snakemake

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
OUT=/DATA/casadei7/nanopore/WGS/RESULTS

python "$BASE/scripts/audit_saved_lrs_outputs.py" \
    --root "$OUT" --sample 11446 --check \
    --report "$BASE/saved_lrs_inventory.tsv"
```

This lists missing required/optional files and exits nonzero if any essential
input is missing. `--sample` can be repeated. Omit it to audit all numeric
patient subfolders. You may move the report outside the repository if it
contains sensitive path information.

## Patient phenotype setup

You still need the secure
`/DATA/casadei7/private_clinical/patient_profiles.tsv` with an entry for
each sample ID, containing only clinically observed `hpo_present` terms.
The independent custom HON ranking retains its core ocular reference
terms. For the primary Genomiser Phenopacket, this wrapper currently sets
`hpo_include_shared_hon_seeds: false`, so Genomiser uses only each patient's
documented HPO findings. Setting it to `true` adds the six core HON reference
terms as computational scoring assumptions; they must not be described as
observed symptoms. A blank clinical field is unknown, not a negative finding.

The installed Exomiser 15.1.0 / 2512 data directories, ReMM v0.4
GRCh38 and tabix index, `config/genomiser/genome_analysis.yml`, HPO
seeds, Monarch nodes/edges, MitoCarta resources, gene BED/panel and VEP
cache referenced in the original LRS config must exist.

## Execute

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/lrs_downstream_reanalysis

# Confirm the selected patients and rule graph
snakemake -s Snakefile --list-rules
snakemake -s Snakefile --cores 8 --dry-run

# Recompute normal saved-data downstream targets, including Genomiser
snakemake -s Snakefile --use-conda --cores 8 \
    --resources mem_mb=45000 --printshellcmds

# Optionally regenerate plots that depend only on saved evidence
snakemake -s Snakefile saved_evidence_plots \
    --use-conda --cores 8 --resources mem_mb=45000
```

The normal target writes the existing final clinical/SV/gene/mitochondrial
tables and filtering outputs beneath each sample folder, and adds
`<sample>/genomiser/` results. To force re-analysis of specific complete
downstream stages, use explicit targets and `--forcerun <rule>`;
do not use `--forceall` indiscriminately. Preserve backup copies of
important outputs before overwriting/recomputing them.

The retained source results must be internally consistent and use the
same GRCh38 coordinate build. A successful command-line test is not
clinical validation; Genomiser support for Jasmine symbolic SV and BND
representation must be checked in the first real sample.

## Cannot be rerun without reads

Fresh 10 kb mosdepth, modkit methylation calling, read-haplotype XCI,
read-based breakpoint figures, fresh SV callers, and read-level
phasing require BAMs or other original read data. This workflow reuses
existing derived outputs wherever possible; it excludes those read-level
operations. Complete `final_thesis_analysis` / `all_thesis_plots`
targets in the original LRS workflow may still require BAMs and must
not be used in the BAM-free mode.

## Optional CADD sensitivity run and unified clinical table

The original Genomiser `config/genomiser/genome_analysis.yml` and the original
clinical ranked TSV are **unchanged**. A second analysis YAML,
`config/genomiser/genome_analysis_cadd.yml`, differs only by adding
`CADD` to `pathogenicitySources`. The optional rule
`genomiser_genome_cadd` executes it on exactly the same VCFs and phenotype
packet as the baseline, but writes to `<sample>/genomiser_cadd/` instead of
`<sample>/genomiser/`. The merger checks sample ID, input paths,
Phenopacket checksum, and different YAML hashes before comparing results.

**Check actual score resources before enabling.** A line such as
`cadd.version=1.7` and **commented-out**
`#exomiser.hg38.cadd-snv-path` /
`#exomiser.hg38.cadd-in-del-path` in the installed
`application.properties` **does not show that either data file is installed**.
CADD requires a GRCh38 SNV file, a compatible GRCh38 indel file and their
matching `.tbi` indexes; indel filenames vary by release.

To locate files on the server, start with:

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
CLI="$BASE/reference/AnnotSV_annotations/jar/exomiser-cli-15.1.0"
find "$CLI/data" /DATA/Reference -type f \
  \( -name 'whole_genome_SNVs.tsv.gz' -o -name 'InDels.tsv.gz' \
     -o -name '*indel*.tsv.gz' \) 2>/dev/null | head -40
```

Then set **only this wrapper's** `config.yaml` to:

```yaml
genomiser_cadd_enabled: true
genomiser_cadd_analysis_yaml: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/config/genomiser/genome_analysis_cadd.yml"
genomiser_cadd_snv_hg38: "/actual/path/to/GRCh38/whole_genome_SNVs.tsv.gz"
genomiser_cadd_indel_hg38: "/actual/path/to/GRCh38/indel-scores.tsv.gz"
```

Replace the **example** score paths with the verified existing absolute
paths; do not type them literally. No changes to the Exomiser
`application.properties` file are required because the CADD
comparison runner supplies those paths through Java system properties.
The baseline Genomiser run is unaffected.

Use the established Snakemake command format from this directory:

```bash
snakemake -s Snakefile --configfile config.yaml --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 --dry-run -p cadd_comparison_all

snakemake -s Snakefile --configfile config.yaml --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 --resources mem_mb=45000 --rerun-incomplete -p cadd_comparison_all
```

When `genomiser_cadd_enabled: false`, the original default `all` target
still runs normally. Enable CADD **only after** verifying resources;
`cadd_comparison_all` is registered when CADD is enabled.

The optional comparison writes one per-patient table, preserving all
original clinical columns and adding `GENOMISER_CADD_*` columns:

```text
RESULTS/<sample>/gene_discovery/final/
  <sample>_clinical_candidates.ranked.tsv
  <sample>_clinical_candidates.genomiser_cadd_comparison.tsv

RESULTS/cohort_analysis/
  clinical_baseline_vs_cadd.tsv
```

The **cohort TSV** is the single combined table across the selected
patients, with `SAMPLE_ID` identifying each row. The comparison has
`GENOMISER_GENE_RANK`, `GENOMISER_CADD_GENE_RANK` and
`GENOMISER_CADD_RANK_IMPROVEMENT` (positive = improved rank with CADD);
also model and SV matching fields for each version. The new CADD columns
are scores/ranks from **Genomiser with CADD included**, not direct CADD
PHRED scores. Do not report them as raw CADD predictions or validated
pathogenicity. Cases absent from a Genomiser run have missing ranks.

The comparison target depends on the original clinical output and can
recompute missing downstream dependencies. Inspect the dry-run carefully,
check the archived-input audit, and preserve a backup of original downstream
results before running. No BAM-based callers or new coverage jobs are included in the archived downstream mode.

## Genomiser Java heap exhaustion from sequence-encoded giant deletions

Exomiser/Jannovar may fail with `java.lang.OutOfMemoryError: Java heap
space` in `NucleotideDeletion.toHGVSString` when caller SV VCF records
encode multi-megabase deletions as full REF sequences (e.g. 48 Mb REF
for a deletion). The Genomiser runner now uses
`scripts/prepare_genomiser_sv_vcf.py` to make a **separate**
`<sample>/genomiser/sv.genomiser.prepared.vcf.gz` (with .tbi index).

For unambiguous single-allele, left-anchored deletion records with
`REF` of at least 1000 bases, it replaces the long REF by its anchor
base and ALT by `<DEL>`, and derives `END` and negative `SVLEN`
from the original REF length. It never modifies the original
`sv/merged/<sample>_merged_SV.vcf.gz`, other SV types, or genotype calls.
Changes are traceable in
`<sample>/genomiser/sv.genomiser.large_deletion_changes.tsv`,
including the original `END`, `SVLEN` and inconsistency flags.
Review flagged discrepancies before any clinical interpretation.

Unexpected unconvertible explicit alleles larger than 1 Mb stop the run
with an actionable message rather than silently passing them to Jannovar.
This does not remove variants or validate their biological plausibility.
The separate custom AnnotSV/clinical SV evidence remains unchanged.
Run `python -m unittest discover -s tests -p test_genomiser_sv_preparation.py`
from the repository root for the small synthetic regression tests.

The separate unrelated breakend `INFO/END < POS` header/interval warning
and conflicting `AF/GL/PL` header definitions are *not fixed* by this
change. They require separate caller-VCF quality review. If a subsequent
Genomiser run still exhausts memory, inspect the new Genomiser Java log;
do not assume SV preparation guarantees success.

## Regenerate all available thesis plots from saved evidence

The `all_saved_thesis_plots` target requests **all available** saved-evidence
plots for every selected patient: Python candidate loci (including detailed
candidate cards), MitoCarta/mitochondrial context, SV and gene overviews,
R phenotype/evidence/population/large-SV figures, optional saved-QC figures,
filtering plots, and cohort comparison only when that module is enabled.
The runner is intentionally **not** called with `--focused`. It cannot
rebuild missing BAM-derived evidence, read-level alignments or XCI.

For a **first run** or normal incremental downstream update:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/lrs_downstream_reanalysis

snakemake -s Snakefile --configfile config.yaml --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs \
  --cores 32 --resources mem_mb=45000 --rerun-incomplete -p \
  all all_saved_thesis_plots
```

**To regenerate every available plot on every invocation, even if it is
already up to date, use the included wrapper**:

```bash
bash run_all_plots.sh
```

This wrapper runs `all all_saved_thesis_plots` with `--forcerun` for the
saved-evidence plotting rules only (Python launcher, downstream filtering,
R phenotype/evidence/population/depth and available R QC). It does **not**
force Genomiser, variant calling, reference generation, or candidate
ranking. Snakemake still recomputes missing/outdated downstream inputs.
The repository is on a `noexec` filesystem on this server, so call the
wrapper with `bash` rather than `./run_all_plots.sh`.
Use `SNAKEMAKE_CONDA_PREFIX`, `SNAKEMAKE_CORES`, and `SNAKEMAKE_MEM_MB`
to override the wrapper's Conda prefix, CPU and memory settings.

Python output figures appear in
`<sample>/plots_reanalysis/` rather than mixing with older
`plots/` results. The figure launcher is independent of removed legacy read-alignment figures. The historical Exomiser-comparison plot is not produced.
Instead `10_integrated_evidence/<sample>_genomiser_rank_comparison.*`
shows descriptive top-N agreement from the **current Genomiser**
comparison tables.

Review `<sample>/plots_reanalysis/plot_run_manifest.tsv` for every
plot script. `COMMAND_COMPLETED` means the script finished, not that
its figure is diagnostically validated; `SKIPPED_MISSING_INPUT` means
that a supporting archived file was absent and must not be silently
interpreted as evidence that the biological phenomenon is absent.
Optional plots require genuine saved outputs: caller support, mosdepth
summaries, needLR, phased VCFs, bedMethyl, genome BED/FAI, etc.
The standard R plots remain in `plots_r/` and cohort comparisons in
`cohort_plots/`.

Publication selection: prioritize SV-type/cohort context, candidate
event evidence, locus breakpoint/depth context, HON phenotype
relevance, and comparison of the two ranking methods. Use repeat,
MEI, phased/methylation and mitochondrial pathway plots when relevant
to an individual candidate, and place broad exploratory figures in
the supplement. Avoid calling multi-caller agreement or phasing
independent experimental validation. Noncoding, regulatory and
inversion-spanned gene context are hypotheses until independently
supported.

## Full available figures, no cohort recurrence, and sensitivity checks

For incremental batches of unresolved patients, `cohort_analysis_enabled: false`
is the default. No SV recurrence is calculated or used by the final ranking;
its status is `NOT_EVALUATED`, **not** absence of recurrence. Shared HON HPO
reference construction is not a cohort statistic and is now stored under
`reference_checks/hon_hpo_reference.tsv`. When the full, comparable cohort
is ready, deliberately switch `cohort_analysis_enabled: true`, confirm
sample selection and rerun `all` using ordinary Snakemake provenance tracking.

`all_saved_thesis_plots` now requests the complete set of currently
supported saved-evidence figures, including:
- SV landscape, candidate evidence matrix, caller agreement and population context
- Mechanism-aware large-event overview and detailed candidate locus/cards
- Phenotype heatmaps, gene prioritization and Genomiser rank comparison
- MitoCarta / nuclear mitochondrial gene and mtDNA context figures when annotated
- Measured coverage, gene/exon and matched caller breakpoint tracks
- R phenotype, evidence, population, saved large-SV depth, QC (when available),
  and filtering plots

Each plot has its own input requirements. Optional figures may be skipped,
and empty folders may still be created; inspect the manifest and log. The
`14_mitochondrial_context/` figures need MitoCarta annotation columns.
All available outputs can be extensive for `candidate_locus_top_n: 50`.

The new compact locus figures and source manifests appear at
`<sample>/plots_reanalysis/12_breakpoint_tracks/`. This panel is
based on archived evidence: it uses the archived 10-kb mosdepth windows and existing
caller VCF interval coordinates. Small events below the depth-bin size
may be unresolved. A dashed marker for a caller appears only if an actual
same-type, same-chromosome record matches both breakpoints; missing
caller VCFs imply unknown support, not a negative result. Do not equate
such a track with split-read/assembly or experimental validation.

The original `gene_bed` contains whole-gene spans only. For actual exon
rectangles set `candidate_plot_gene_gtf` in the inherited LRS
`config_lrs.yaml` to an existing GRCh38 GENCODE GTF path. With
`null`, the plot correctly shows gene spans but **does not invent exons**.

The standard `all` target also writes the non-destructive
`ranking_stability/` table and plot by comparing the unmodified final
event ordering under alternative population-frequency and multi-caller
eligibility filters. It is not a truth-based diagnostic benchmark and it
does not rerun or rescore the ranking. Unknown population AF is explicitly
retained for review, not treated as AF zero. Configure
`ranking_stability_enabled: false` to skip this optional sensitivity check.

The shared `benchmarking_review_all` target writes a concise per-patient
`benchmark/<sample>.benchmark_review.tsv` and JSON manifest. It links
the existing descriptive method agreement, eligibility-filter ranking
stability and (only when independently curated) confirmed-truth results.
There is no new pathogenicity score or patient RNA-expression analysis.

Storage cleanup is performed separately and is never part of any Snakemake target. Preserve source evidence for reproducible replotting.
