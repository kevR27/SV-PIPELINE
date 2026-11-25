# LRS downstream reanalysis from saved outputs (without BAM)

This directory runs the **existing LRS downstream interpretation** from
previously generated files. It reuses the rules in
`../lrs/Snakefile_LRS_postprocess`, instead of maintaining a separate,
possibly divergent copy of the scoring algorithms. It **does not rerun**
basecalling, Clair3, SV callers, Jasmine, LongPhase, mosdepth, modkit or
Samplot. Genomiser is included and mandatory.

**Different from `../genomiser_vcf_only/`:** this workflow rebuilds with retained BAM-derived evidence
 the full saved-data SV gene/event prioritization, multimodal integration,
mitochondrial ranking, cohort recurrence, clinical candidate TSVs, downstream
filtering and selected R/summary plots. It therefore needs more of the
previous run's BAM-derived TSV/VCF/BED output files.

## Automatic sample discovery

Set the common root only once:

```yaml
input_root: "/DATA/casadei7/nanopore/WGS/RESULTS"
sample_ids: auto
```

Numeric subfolders such as `11446/`, `21930/`, `25428/` are detected
automatically. Original filenames are derived from the sample directory
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
The six shared HON HPO terms are added for computational scoring when
clinically documented optic atrophy `HP:0000648` is present. Reference
terms do not become claimed patient-observed symptoms.

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
read-based SV support plots (Samplot), fresh SV callers, and read-level
phasing require BAMs or other original read data. This workflow reuses
existing derived outputs wherever possible; it excludes those read-level
operations. Complete `final_thesis_analysis` / `all_thesis_plots`
targets in the original LRS workflow may still require BAMs and must
not be used in the BAM-free mode.
