# Patient-specific optic neuropathy analysis using Genomiser

## Status

This is an **opt-in implementation**. The Snakemake workflow structure has
been updated, but it has **not** been validated on the laboratory server,
real ONT/Illumina VCFs, or the local ReMM database. Keep
`genomiser_enabled: false` until the checks below succeed.

**Genomiser** is the `--preset genome` analysis in **Exomiser CLI 15.1.0**.
It analyses coding and noncoding variants jointly with the phenotype. It is
not a separate program named "Genomizer". It requires a matching GRCh38 ReMM
score file and tabix index in addition to the existing 2512 genome and
phenotype resources. See:
- https://exomiser.readthedocs.io/en/stable/input_files_and_options.html
- https://exomiser.readthedocs.io/en/stable/sv_prioritisation.html
- https://exomiser.readthedocs.io/en/stable/installation.html
- https://remm.bihealth.org/download

The original phenotype-only Exomiser comparison is no longer a mandatory
`rule all` target. Existing Exomiser scripts are retained for historical
reproducibility; they do not run by default.

## One private patient table for both LRS and SRS

Copy `reference/patient_profiles.template.tsv` to a secure location,
for example `/DATA/casadei7/private_clinical/patient_profiles.tsv`.
**Do not put patient clinical records in this public GitHub repository.**
Keep one tab-separated row per sample, with the following columns:

- `sample_id`: exactly the ID in `samples:` in the corresponding YAML.
- `hpo_present`: semicolon-separated, *documented* HPO IDs.
- `hpo_absent`: optional, explicitly examined and absent terms.
- `onset_min_years`, `onset_max_years`: numeric age range; leave blank if unknown.
- `onset_status`: `confirmed`, `estimated`, `unknown`.
- `sex`: `MALE`, `FEMALE`, `UNKNOWN_SEX` if known (important for sex-linked MOI).
- `diagnostic_hypothesis` and `hypothesis_status`: reviewer annotations, not
  confirmed genetic diagnoses or scoring features.

The six `CORE_OCULAR_HON` terms in `reference/hon_hpo_seed_terms.tsv` are
**shared reference ontology seeds only**. Do **not** assign all six to each
patient unless the patient actually exhibits each finding.

Additional neurological/syndromic findings should remain in a given
patient's `hpo_present` when clinically documented, even though generic
mitochondrial/syndromic seeds were removed from the cohort template.

The importer generates:
```
RESULTS/<sample>/clinical/<sample>.patient_hpo.tsv
RESULTS/<sample>/clinical/<sample>.phenopacket.json
RESULTS/<sample>/clinical/<sample>.clinical_context.json
```
Both your final custom phenotype comparison and Genomiser use the same
**present** HPO set. Absent HPO, age-at-onset estimates and clinical diagnostic
hypotheses are retained for independent manual review, *not* made into
numerical priors for either algorithm. They therefore cannot leak the
expected disease into the benchmark. Ontology IDs are syntax-validated;
verify the IDs and their meanings against the current HPO release before
clinical use.

## Configure

In BOTH `snakemake_pipelines/lrs/config_lrs.yaml` and
`snakemake_pipelines/srs_wgs_pipeline/config_srs_wgs.yaml`, set:

```yaml
patient_profiles_tsv: "/DATA/casadei7/private_clinical/patient_profiles.tsv"
genomiser_enabled: true
genomiser_remm_hg38: "/DATA/casadei7/reference/ReMM.v0.4.hg38.tsv.gz"
genomiser_memory_gb: 32
genomiser_benchmark_truth_tsv: null
```

Verify that `genomiser_cli_jar`, `genomiser_data_dir`, the existing
`2512_hg38` and `2512_phenotype` directories, the ReMM file and
`ReMM.v0.4.hg38.tsv.gz.tbi` all exist. Keep a reserve of RAM above the Java
heap, avoid simultaneous 32GB jobs if memory is limited, and check free disk
space before combining genome-wide VCFs.

To validate profiles without any sequencing run:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3
python scripts/prepare_patient_profiles.py \
  --profiles /DATA/casadei7/private_clinical/patient_profiles.tsv \
  --sample 21930 \
  --hpo-seeds reference/hon_hpo_seed_terms.tsv \
  --phenopacket /tmp/21930.phenopacket.json \
  --hpo-output /tmp/21930.patient_hpo.tsv \
  --clinical-output /tmp/21930.clinical_context.json
```

Then, from the LRS folder:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/lrs
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --dry-run
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --resources mem_mb=65000
snakemake -s Snakefile_LRS_postprocess --use-conda --cores 32 --dry-run
snakemake -s Snakefile_LRS_postprocess --use-conda --cores 32
```

For SRS:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/srs_wgs_pipeline
snakemake -s Snakefile_SRS_WGS --use-conda --cores 32 --dry-run
snakemake -s Snakefile_SRS_WGS --use-conda --cores 32 --resources mem_mb=65000
```

Adapt the total memory resource to the actual server. Keep all existing
calling, filtering and annotation outputs; the new genome analysis is
additive. Running `--dry-run` is not sufficient to prove the Java analysis
works. Start by inspecting one positive control, its combined VCF, and
its Genomiser report.

## What is compared

Per sample, Genomiser receives the original genome-wide SNV/indel VCF
(Clair3 on ONT; DeepVariant on Illumina) plus a separate SV VCF (Jasmine
on ONT; integrated Manta/Delly/CNV on Illumina). The script combines these
into a single sample-matched, sorted and indexed VCF while preserving the
original files. Genomiser writes its own `genes.tsv`, `variants.tsv`, and
HTML report under `RESULTS/<sample>/genomiser/`.

LRS additional evidence such as methylation, long-read phasing, tandem
repeats and MEIs and SRS additional evidence such as repeat expansion,
MELT and independent read depth do **not** become Genomiser inputs merely
because they appear in the custom integrated TSV. Explicitly check variant
representation of symbolic DEL/DUP/INV/INS and breakends: Genomiser's SV
support has been tested on selected SRS callers, and an ONT Jasmine VCF
requires compatibility validation before trusting absent results.

## Independent benchmarking (only when clinical truth exists)

Make a **private**, pre-curated file following
`reference/benchmark_truth.template.tsv`; it must contain `sample_id`,
`gene`, `truth_status=CONFIRMED`, and `variant_class`. For a breakpoint-
validated DEL/DUP/INS/INV, also include `chrom`, `pos`, `end`,
`svtype`. Otherwise use `variant_class=GENE_ONLY`. Put only independently
validated truth in this file; a clinical suspicion or algorithmic call is
not ground truth.

Set `genomiser_benchmark_truth_tsv` to the private file in both YAML
configs. In LRS, the benchmark runs in
`Snakefile_LRS_postprocess`; in SRS, it runs in
`Snakefile_SRS_WGS`. Outputs are per-sample:
```
RESULTS/<sample>/benchmark/<sample>.genomiser_vs_custom.tsv
RESULTS/<sample>/benchmark/<sample>.genomiser_vs_custom.json
```
Gene-level top-K recall is the proportion of independently solved truth
entries with the causal gene in the top K. The SV-level comparison uses
a provisional 500-bp positional tolerance for supported symbolic types
and must be manually checked. The Genomiser variant TSV's `#RANK` is
associated with its gene/MOI score; it must **not** be presented as an
independent per-SV rank.

For a valid method comparison, stratify by LRS/SRS, variant class,
clinical completeness, and candidate callability. Separate *variant
detection* (truth allele missing from the input VCF) from *ranking failure*
(allele present but not prioritized). Without independently solved controls,
rank correlation and discordance can be reported, but not sensitivity or
diagnostic superiority. Preserve the original ranking; no weights are
trained/changed using the benchmark truth.
