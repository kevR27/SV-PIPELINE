# Patient-specific optic neuropathy analysis using Genomiser

## Status

**Genomiser is now mandatory** for every patient in the normal LRS
`Snakefile_LRS_update`, LRS `Snakefile_LRS_postprocess`, and SRS
`Snakefile_SRS_WGS` analyses. Their default `all` targets require its
per-sample manifest; the separate thesis plot/analysis targets also require
it. Setting `genomiser_enabled: false` raises an explicit configuration error.
The LRS and SRS final candidate-ranking rules also depend on the
Genomiser manifest, so directly requesting an updated final ranked TSV
cannot bypass the required genome analysis. Genomiser does not alter
the custom ranking algorithm; this is an execution dependency only.

**This does not mean the laboratory server has been validated.** Before
running either workflow, install the GRCh38 ReMM database/index and populate
the secured patient profile table. Missing prerequisites now stop the pipeline;
they no longer allow analysis to finish without Genomiser.

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


## Editable analysis configuration (required)

The code now uses an **explicit YAML** passed to Exomiser 15.1.0 with
`--analysis`, **not** `--preset genome`. The repository includes a copy
of the upstream 15.1.0 **genome** preset:

```text
config/genomiser/genome_analysis.yml
```

The baseline has not been altered biologically; it includes coding,
intronic, intergenic and regulatory variant evaluation using REVEL,
MVP, AlphaMissense, ReMM and SpliceAI as available in the installed
databases. It retains the upstream `PASS_ONLY` mode, phenotype-first
priority filter and genome-wide regulatory and pathogenicity filters.
This is a starting baseline, not an assertion that every input variant
will survive every downstream filter.

After `git pull`, edit directly on the server:

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
nano "$BASE/config/genomiser/genome_analysis.yml"
```

Both LRS and SRS YAML configs have the required entry:

```yaml
genomiser_analysis_yaml: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/config/genomiser/genome_analysis.yml"
```

The runner loads this file with PyYAML, checks that it is a genuine
genome-wide analysis with ReMM and population-frequency assessment,
and refuses gene-panel/interval restrictions or a variantEffectFilter
that removes intronic, regulatory or intergenic variant classes.
Other analysis options can be changed in the YAML. **Every edit to the
analysis YAML is a declared Snakemake input change**, so the
per-patient Genomiser analysis will be rerun automatically on the
next requested downstream target. Each result manifest records the
exact analysis YAML SHA-256 and its selected steps and sources.

Change one parameter at a time when benchmarking against solved
controls. `analysisMode: FULL` can be used to retain failed variants
for troubleshooting, but requires substantially more RAM/disk than
`PASS_ONLY`. The runner keeps gene and variant TSV output enabled and
writes `output-options.yml` automatically: editing the per-sample
generated copy will *not* change a later run.

For source and documentation:
- Baseline: https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-cli/src/main/resources/examples/preset-genome-analysis.yml
- Configuration: https://exomiser.readthedocs.io/en/stable/advanced_analysis.html
- CLI `--analysis`: https://exomiser.readthedocs.io/en/stable/input_files_and_options.html

## Reanalysis starting from the original VCFs

**Recommended:** keep the original calling results, and rerun the
variant-interpretation stage against those same files. Do not delete
the VCFs or use `--forceall` on the complete pipeline unless you
intend to repeat expensive upstream calling.

Per patient, the Genomiser runner needs:

1. Genome-wide small-variant `VCF.gz`: LRS Clair3
   `RESULTS/<sample>/snp_clair3/<sample>.vcf.gz` or SRS
   DeepVariant `RESULTS/<sample>/snp_deepvariant/<sample>.vcf.gz`.
2. Genome-wide SV `VCF.gz`: LRS Jasmine-merged
   `RESULTS/<sample>/sv/merged/<sample>_merged_SV.vcf.gz`
   or SRS `RESULTS/<sample>/sv/integrated/<sample>_integrated_SV.vcf.gz`.
3. The private `patient_profiles.tsv` containing the patient's
   positive HPO terms; the pipeline builds a Phenopacket automatically.
4. `config/genomiser/genome_analysis.yml` (editable analysis filters).
5. Exomiser CLI 15.1.0 JAR with the GRCh38/phenotype 2512 databases
   and correctly configured `application.properties`; indexed
   ReMM v0.4 hg38 data.
6. `envs/genomiser.yaml` (Snakemake creates or reuses the environment).

Do not fabricate technical-control phenotypes. Note that the more extensive
custom SV ranking can require its **existing** AnnotSV, read-depth, VEP,
phasing and orthogonal-evidence tables if you intend to rebuild it too.
Keep or regenerate those downstream dependencies separately.

Once source VCFs exist, check Genomiser alone for sample 21930:

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
OUT=/DATA/casadei7/nanopore/WGS/RESULTS
SAMPLE=21930

cd "$BASE/snakemake_pipelines/lrs"
snakemake -s Snakefile_LRS_update -n --cores 32 \
  "$OUT/$SAMPLE/genomiser/$SAMPLE.genomiser.manifest.json"
```

If the dry-run job list does **not** include genome variant calling,
run the same command without `-n` and with `--use-conda`. Running
the LRS postprocess final ranked TSV target afterwards ensures that
the custom ranking is rebuilt when its inputs have changed.
For SRS, use `Snakefile_SRS_WGS` and the output directory and
sample ID configured in `config_srs_wgs.yaml`, not the LRS paths.

To deliberately repeat Genomiser with **unchanged** YAML and VCFs,
use `--forcerun genomiser_genome` on the requested sample output.
This forces that analysis only, not the original callers.

## Required installation / resource verification on the Bologna server

Because the 15.1.0 Exomiser CLI and its 2512 HG38/phenotype databases are
already configured for the existing phenotype-only comparison, **do not
download duplicate genome databases unless the preflight identifies missing
files**. Retain 15.1.0 and 2512 for reproducibility; a newer Exomiser release
requires a separate versioned validation.

1. Confirm the installed Java CLI and data:

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
CLI="$BASE/reference/AnnotSV_annotations/jar/exomiser-cli-15.1.0"
test -s "$CLI/exomiser-cli-15.1.0.jar"
find "$CLI/data" -maxdepth 2 -type f -name '*mv.db' | head -20
df -h "$BASE"
```

2. Download ReMM v0.4 **GRCh38** (approximately 15 GB compressed) and the
index from the official ReMM page (https://remm.bihealth.org/download).
Use the direct official URLs below; `wget -c` supports continuing an
interrupted download:

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
REMM_DIR="$BASE/reference/remm"
mkdir -p "$REMM_DIR"
cd "$REMM_DIR"
wget -c https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz
wget -c https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.tsv.gz.tbi
wget -c https://kircherlab.bihealth.org/download/ReMM/ReMM.v0.4.hg38.md5
md5sum -c ReMM.v0.4.hg38.md5
ls -lh ReMM.v0.4.hg38.tsv.gz*
```

If any checksum fails, stop and redownload the corrupted file; do not
create a new tabix index as a substitute for a missing official file.

3. Install the small execution environment using the workflow's YAML:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3
conda env create -f envs/genomiser.yaml
conda run -n genomiser java -version
conda run -n genomiser bcftools --version
```

If `genomiser` already exists, use
`conda env update -n genomiser -f envs/genomiser.yaml --prune`
instead of creating it. Snakemake `--use-conda` can create its own
environment from the same file, so a manual environment is for standalone
preflight or debugging.

4. Edit the **two** YAML configs:
```yaml
patient_profiles_tsv: "/DATA/casadei7/private_clinical/patient_profiles.tsv"
genomiser_enabled: true
genomiser_analysis_yaml: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/config/genomiser/genome_analysis.yml"
genomiser_remm_hg38: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/remm/ReMM.v0.4.hg38.tsv.gz"
genomiser_memory_gb: 32
genomiser_benchmark_truth_tsv: null
```

Check `genomiser_cli_jar` and `genomiser_data_dir` point to the existing
installation. The runner passes the 2512 data versions and ReMM path via
Java system properties; the official `application.properties` in the
installation directory is loaded when present.

5. Add the clinical table, ensuring **HP:0000648** is listed in
`hpo_present` for every OA patient. Include *only confirmed observed*
additional findings. In particular, missing or unexamined features are
not negative HPO observations:

```bash
mkdir -p /DATA/casadei7/private_clinical
chmod 700 /DATA/casadei7/private_clinical
cp "$BASE/reference/patient_profiles.template.tsv" \
   /DATA/casadei7/private_clinical/patient_profiles.tsv
nano /DATA/casadei7/private_clinical/patient_profiles.tsv
chmod 600 /DATA/casadei7/private_clinical/patient_profiles.tsv
```

Replace the example row with actual de-identified patient sample IDs
matching the YAML `samples` keys.

## Rerun ONLY phenotype ranking and Genomiser

**Do not invoke the default `all` target simply to update patient HPO
or the Genomiser comparison**: Snakemake may rebuild missing dependencies
from expensive callers. First run `--dry-run` for exact targets.

For LRS example patient `21930`:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/lrs
OUT=/DATA/casadei7/nanopore/WGS/RESULTS
snakemake -s Snakefile_LRS_update --cores 32 --dry-run \
  "$OUT/21930/genomiser/21930.genomiser.manifest.json"
snakemake -s Snakefile_LRS_update --use-conda --cores 32 \
  --resources mem_mb=65000 \
  "$OUT/21930/genomiser/21930.genomiser.manifest.json"
snakemake -s Snakefile_LRS_postprocess --cores 32 --dry-run \
  "$OUT/21930/gene_discovery/final/21930_gene_candidates.ranked.tsv"
snakemake -s Snakefile_LRS_postprocess --use-conda --cores 32 \
  "$OUT/21930/gene_discovery/final/21930_gene_candidates.ranked.tsv"
```

For SRS, use the actual `path` and `sample` in
`config_srs_wgs.yaml`; they can differ from LRS:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/srs_wgs_pipeline
SAMPLE=HG002
OUT=/DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/srs_wgs_pipeline/srs_CNTRL_results
snakemake -s Snakefile_SRS_WGS --cores 32 --dry-run \
  "$OUT/$SAMPLE/genomiser/$SAMPLE.genomiser.manifest.json" \
  "$OUT/$SAMPLE/gene_discovery/final/${SAMPLE}_gene_candidates.ranked.tsv"
snakemake -s Snakefile_SRS_WGS --use-conda --cores 32 \
  --resources mem_mb=65000 \
  "$OUT/$SAMPLE/genomiser/$SAMPLE.genomiser.manifest.json" \
  "$OUT/$SAMPLE/gene_discovery/final/${SAMPLE}_gene_candidates.ranked.tsv"
```

HG002 is a technical example, not a hereditary optic neuropathy patient,
and must not be supplied with fabricated OA HPO observations. Replace
`SAMPLE` and `OUT` with a genuine OA sample and its configured SRS result
directory when running patient-specific clinical prioritization.

When a dry run lists expensive calling jobs, check whether source outputs
were deleted or paths changed. Restore retained source VCFs where possible
rather than recomputing the entire pipeline.

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

## Required configuration

The shipped YAML files already specify the secured patient table and ReMM
paths, but those files must exist on the server before execution.
In BOTH `snakemake_pipelines/lrs/config_lrs.yaml` and
`snakemake_pipelines/srs_wgs_pipeline/config_srs_wgs.yaml`, set:

```yaml
patient_profiles_tsv: "/DATA/casadei7/private_clinical/patient_profiles.tsv"
genomiser_enabled: true
genomiser_remm_hg38: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/remm/ReMM.v0.4.hg38.tsv.gz"
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

For a normal full re-analysis, `rule all` always requires Genomiser.
To update only the phenotype ranking and Genomiser, use the targeted
commands above: sequencing and variant calling do not need to be repeated
when their outputs are intact.

All configured OA patients must have a genuine clinical row with
`HP:0000648` in `hpo_present`; never invent an optic atrophy
phenotype for a technical control such as HG002. The default SRS HG002
sample is a test and cannot be analysed as an OA patient unless you replace
it with an appropriately phenotyped patient or supply legitimate clinical
observations for the control itself.

**Benchmarking remains separately conditional on independently confirmed
diagnostic truth.** The Genomiser stage itself is never conditional on truth.

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
Gene-level top-K recall is evaluated using a **reconstructed global gene order** for the custom ranking (independent of panel status) and the global Genomiser gene order. For independently solved SVs, the script also reports whether the corresponding allele appears in each method's output among the top-K genes. A provisional 500-bp positional tolerance is used for simple symbolic SV matches, which must be manually checked. The Genomiser variant TSV's `#RANK` corresponds to its gene/MOI model, **not** a separately ranked SV. Unsolved samples are skipped by the optional benchmark target.

For a valid method comparison, stratify by LRS/SRS, variant class,
clinical completeness, and candidate callability. Separate *variant
detection* (truth allele missing from the input VCF) from *ranking failure*
(allele present but not prioritized). Without independently solved controls,
rank correlation and discordance can be reported, but not sensitivity or
diagnostic superiority. Preserve the original ranking; no weights are
trained/changed using the benchmark truth.
