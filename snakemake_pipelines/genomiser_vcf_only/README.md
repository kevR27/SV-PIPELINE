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

Example for `21930`:

```text
/DATA/casadei7/nanopore/WGS/RESULTS/21930/
  snp_clair3/
    21930.vcf.gz
    21930.vcf.gz.tbi
  sv/
    merged/
      21930_merged_SV.vcf.gz
      21930_merged_SV.vcf.gz.tbi
```

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

Replace the example sample with the actual identifiers in
`vcf_samples`. One row per patient; include documented `hpo_present`
and, if available, explicitly assessed `hpo_absent`. The pipeline
automatically adds the shared HON seed terms to computational scoring
for patients with documented optic atrophy (`HP:0000648`).
Do not fabricate phenotypes or copy private clinical data into GitHub.

## Run

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/genomiser_vcf_only

nano config.yaml

snakemake -s Snakefile --cores 8 --dry-run
snakemake -s Snakefile --use-conda --cores 8 \
  --resources mem_mb=45000 --printshellcmds
```

The run creates `<path>/<sample>/clinical/` and
`<path>/<sample>/genomiser/`. It does not generate the original
LRS/SRS custom-integrated clinical ranking or perform new variant calling.

For more detail, see `../../docs/GENOMISER_VCF_ONLY.md`.
