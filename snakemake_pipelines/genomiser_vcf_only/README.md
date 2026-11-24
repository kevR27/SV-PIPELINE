# VCF-only Genomiser

Standalone, BAM-free Genomiser analysis for existing LRS or SRS genome-wide
small-variant and merged SV VCFs. It is separate from the LRS/SRS calling
workflows and cannot schedule their BAM-dependent rules.

## Files in this directory

- `Snakefile`: orchestration (patient profile + Genomiser only).
- `config.yaml`: input VCF paths, output path and Exomiser/ReMM settings.
- `README.md`: these instructions.

Shared code (from the same Git repository) is reused:
`../../scripts/prepare_patient_profiles.py`,
`../../scripts/run_genomiser.py`,
`../common/patient_profiles.smk`,
`../common/genomiser.smk`, and `../../envs/genomiser.yaml`.
Do not copy the entire LRS Snakemake directory or any BAMs here.

## Upload only the required variant files for each sample

Example for the clinical HON patient `14996` (the other active
patient is `13048`):

```text
/DATA/casadei7/nanopore/WGS/RESULTS/14996/
  snp_clair3/
    14996.vcf.gz
    14996.vcf.gz.tbi
  sv/
    merged/
      14996_merged_SV.vcf.gz
      14996_merged_SV.vcf.gz.tbi
```

Validation controls `21930` and `25428` are **not** in this
workflow's `vcf_samples`. The original LRS/SRS config files were
restored and are not modified by this separate analysis.

The two `.vcf.gz` files are mandatory. Matching indexes are
recommended (or use `.csi` when appropriate). Do not supply panel-only
or heavily filtered review VCFs as a substitute for genome-wide calls.
If the files already exist on the server, do not reupload them.
Set `vcf_samples` in `config.yaml` to the actual locations, so
the original result folders need not be copied or renamed.

The original BAM/BAI, separate Sniffles/cuteSV/DELLY VCFs, methylation,
depth, LongPhase, AnnotSV and plot folders are not needed to run Genomiser.
These are needed for other downstream analyses, not this workflow.

## Create the private patient table

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
mkdir -p /DATA/casadei7/private_clinical
chmod 700 /DATA/casadei7/private_clinical
cp -n "$BASE/reference/patient_profiles.template.tsv" \
  /DATA/casadei7/private_clinical/patient_profiles.tsv
chmod 600 /DATA/casadei7/private_clinical/patient_profiles.tsv
nano /DATA/casadei7/private_clinical/patient_profiles.tsv
```

Use the corrected clinical `patient_profiles.tsv`, with `HP:0000648`
confirmed for all 10 HON patients. Missing fields are unknown, not
negative observations. Only explicitly documented absent findings belong
in `hpo_absent`. This VCF-only config intentionally sets
`hpo_include_shared_hon_seeds: true` to retain the six shared
computational HON reference terms alongside confirmed patient findings.
The shared terms are **not** additional clinically confirmed findings.
Do not fabricate phenotypes or copy private clinical data into GitHub.

## Run using the same Snakemake options as the original pipeline

Update `main` on the server, confirm the existing VCF input paths and
private patient table, then run from this **separate** directory.

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/genomiser_vcf_only

snakemake -s Snakefile --configfile config.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 --dry-run -p

# Run only after the dry run schedules the expected patient-profile and
# Genomiser jobs, and only after verifying the VCFs/Exomiser/ReMM resources:
snakemake -s Snakefile --configfile config.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 --resources mem_mb=45000 -p
```

The `path` key points to
`/DATA/casadei7/nanopore/WGS/GENOMISER_VCF_ONLY_RESULTS`
so output files are kept **separate** from
`/DATA/casadei7/nanopore/WGS/RESULTS`.
Only patients `14996` and `13048` are configured for now;
add further clinical patients as their original genome-wide VCFs become available.

The run creates `<path>/<sample>/clinical/` and
`<path>/<sample>/genomiser/`. It does not generate the original
LRS/SRS custom-integrated clinical ranking or perform new variant calling.

For more detail, see `../../docs/GENOMISER_VCF_ONLY.md`.
