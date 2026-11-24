# Running Genomiser from existing VCFs without any BAM

Use the **standalone** `snakemake_pipelines/genomiser_vcf_only/Snakefile` for reanalysis when
the LRS/SRS calling work is already complete and only the genome-wide
SNV/indel VCF and merged SV VCF were retained or restored.

This entry point contains only three rules: `all`, `prepare_patient_profile`,
and `genomiser_genome`. **No BAM, BAI, Clair3, Jasmine, read-depth,
phasing, methylation, needLR, samplot, or caller jobs can be triggered.**
Do not disable rules in `Snakefile_LRS_update`; use this entry point instead.

**Scope:** Genomiser gene/variant TSVs, HTML and manifest, plus patient
Phenopacket/HPO and clinical provenance outputs. This does **not**
reconstruct the full custom SV ranking or read-level evidence. That
requires additional AnnotSV/VEP/caller/clinical context outputs, and some
existing postprocess rules depend on BAMs.

## Directory layout

The standalone Snakemake workflow is kept outside both LRS and SRS,
under `snakemake_pipelines/genomiser_vcf_only/`:

```text
snakemake_pipelines/
  genomiser_vcf_only/
    Snakefile
    config.yaml
    README.md
  common/
    patient_profiles.smk
    genomiser.smk
```

Run it from its own directory. The shared code and Conda environment
remain in the main repository, and existing patient result folders
are unchanged.

## 1. Inputs to restore

For each patient, restore the original **genome-wide** VCFs, not the
panel-only or filtered review VCFs. Example sample `21930`:

```text
/DATA/casadei7/nanopore/WGS/RESULTS/21930/
  snp_clair3/21930.vcf.gz
  snp_clair3/21930.vcf.gz.tbi
  sv/merged/21930_merged_SV.vcf.gz
  sv/merged/21930_merged_SV.vcf.gz.tbi
```

The tabix indexes are strongly recommended; `.csi` indexes are
also valid when appropriate. If the VCFs are properly position-sorted
BGZF files, they can be indexed using `bcftools index -f -t`.
The two VCFs must each have **exactly one genotype sample** and use a
compatible GRCh38 reference and chromosome naming.

**One merged SV VCF alone is not enough for the combined SNV+SV analysis.**
If the small-variant VCF was deleted, restore it; do not invent an
empty small-variant VCF to simulate whole-genome SNV calling.

Required non-VCF resources:
- Existing Exomiser CLI 15.1.0 JAR and extracted `2512_hg38` and
  `2512_phenotype` database directories
- GRCh38 `ReMM.v0.4.hg38.tsv.gz` and `.tbi`
- `config/genomiser/genome_analysis.yml`
- `reference/hon_hpo_seed_terms.tsv`
- Private `/DATA/casadei7/private_clinical/patient_profiles.tsv`
- `envs/genomiser.yaml` plus the other repository code (use `git pull`)

You do not need to upload individual caller BAMs or the separate
Sniffles/cuteSV/DELLY VCFs.

## 2. Create the private patient table

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
CLINICAL=/DATA/casadei7/private_clinical
mkdir -p "$CLINICAL"
chmod 700 "$CLINICAL"

# Only copy the example if a private table does not already exist.
if [ ! -f "$CLINICAL/patient_profiles.tsv" ]; then
  cp "$BASE/reference/patient_profiles.template.tsv" "$CLINICAL/patient_profiles.tsv"
fi

chmod 600 "$CLINICAL/patient_profiles.tsv"
nano "$CLINICAL/patient_profiles.tsv"
```

Use a true tab-separated table with one row per patient and the column
names from `reference/patient_profiles.template.tsv`. The shipped example
`EXAMPLE_001` must be replaced by actual de-identified patient IDs.
Include **only clinically documented findings** in `hpo_present`.
For patients with documented optic atrophy, include `HP:0000648`;
the pipeline augments the Genomiser query with the six shared HON
reference HPO terms and patient-specific terms, excluding documented
absences. The added seeds remain computational reference assumptions,
not clinically verified observations.

## 3. Configure the VCF-only entry point

```bash
cd "$BASE/snakemake_pipelines/genomiser_vcf_only"
nano config.yaml
```

The `vcf_samples` field maps each sample ID to its actual VCF paths.
The `path` field controls the output directory. All other Genomiser
resources have the same defaults as the LRS/SRS workflows; adjust if
installed elsewhere. Edit `vcf_samples` to add patients or SRS VCFs.
There is **no BAM configuration field**.

## 4. Preflight check

```bash
BASE=/DATA/casadei7/tools/SV-PIPELINE-main_v3
OUT=/DATA/casadei7/nanopore/WGS/RESULTS
SAMPLE=21930
bcftools query -l "$OUT/$SAMPLE/snp_clair3/$SAMPLE.vcf.gz"
bcftools query -l "$OUT/$SAMPLE/sv/merged/${SAMPLE}_merged_SV.vcf.gz"

python "$BASE/scripts/prepare_patient_profiles.py" \
  --profiles /DATA/casadei7/private_clinical/patient_profiles.tsv \
  --sample "$SAMPLE" \
  --hpo-seeds "$BASE/reference/hon_hpo_seed_terms.tsv" \
  --phenopacket "/DATA/casadei7/private_clinical/${SAMPLE}.phenopacket.json" \
  --hpo-output "/DATA/casadei7/private_clinical/${SAMPLE}.patient_hpo.tsv" \
  --clinical-output "/DATA/casadei7/private_clinical/${SAMPLE}.clinical_context.json"
```

Verify that the HPO TSV and Phenopacket contain the intended combined
phenotypes before running Genomiser.

## 5. Run without BAMs

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/genomiser_vcf_only

# Verify the job graph: only patient-profile and Genomiser jobs.
snakemake -s Snakefile --cores 8 --dry-run

# Analyse the configured VCF samples; budget one 32-GB Java job at a time.
snakemake -s Snakefile --use-conda \
  --cores 8 --resources mem_mb=45000 --printshellcmds
```

The Genomiser result files are written under
`<path>/<sample>/genomiser/`, including
`<sample>.genomiser.genes.tsv`,
`<sample>.genomiser.variants.tsv` and
`<sample>.genomiser.manifest.json`. The ClinVar/Regulatory database
and merged VCF compatibility must still be validated on the lab server.
For failure details inspect
`<path>/<sample>/logs/genomiser/<sample>.snakemake.log` and
`<path>/<sample>/genomiser/<sample>.genomiser.log`.

If you later need the complete custom final ranking, do not assume these
four VCF files are sufficient: the broader postprocess workflow uses
read-depth, VEP/AnnotSV, orthogonal caller evidence, and other analyses.
