# BAM-free VCF-only Genomiser and custom SV ranking

Use the separate `snakemake_pipelines/genomiser_vcf_only/Snakefile` with its
own `config.yaml`. Both genome-wide VCFs (SNV/indel and merged SV, GRCh38,
single sample each) are required. No BAM callers or read-depth/phasing
producers are registered.

The workflow runs Genomiser and the original shared custom SV ranking scripts:
AnnotSV, VEP, human Monarch gene–HPO evidence, core HON gene ranking,
categorical allele assessment, gnomAD-SV, MitoCarta, patient-aware HPO,
SNV–SV second-allele review and clinical ranking. No read-derived evidence is
fabricated: original per-caller genotype QC, phase, coverage, methylation,
MEI/STR and recurrence are unavailable in a strict two-VCF run. This is not
numerically equivalent to a full BAM-supported analysis.

Configure the same annotation resources, Genomiser YAML, ReMM, private patient
profile TSV and HPO query mode as the standard workflows. This standalone VCF-only workflow is configured with
`hpo_include_shared_hon_seeds: true`, retaining six shared core HON reference
terms **alongside** documented patient HPO findings for computational scoring.
The shared references are not additional observed patient findings; blanks
remain unknown, not absent. The normal LRS/SRS configs currently use `false`.
Do not interpret cross-workflow ranking differences as algorithmic differences
without harmonising these HPO inputs for a matched comparison.

The explicit VCF map can contain LRS or SRS sample paths and a per-sample
`sequencing_type: LRS` or `SRS` (default LRS). Keep the output path distinct
from original BAM-derived results. Configure only the patients being studied.

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/genomiser_vcf_only
snakemake -s Snakefile --configfile config.yaml --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 --dry-run -p
snakemake -s Snakefile --configfile config.yaml --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 \
  --resources mem_mb=45000 --rerun-incomplete -p
```

Results under `<path>/<sample>/` include:

- `genomiser/<sample>.genomiser.genes.tsv`, `*.variants.tsv` and manifest.
- `gene_discovery/final/<sample>_gene_candidates.ranked.tsv`, SV–gene and
  clinical ranked TSVs. The clinical table adds independent Genomiser gene/model
  ranks, native scores when present, and conservative SV coordinate-match status.
- `benchmark/<sample>.ranking_comparison.tsv` and JSON: descriptive global gene
  ranks and top-k overlap for every case, including unresolved cases.
- `benchmark/<sample>.genomiser_vs_custom.tsv` and JSON when an independently
  confirmed truth record is provided through `genomiser_benchmark_truth_tsv`.

Use the header in `reference/benchmark_truth.template.tsv` for a private truth
file. Do not label unresolved hypotheses as confirmed. Top-k overlap is not a
diagnostic accuracy estimate; the two approaches search different candidate
universes. Missing Genomiser results do not prove that a variant was not called.

For a complete reanalysis using saved caller/coverage/phase evidence, use the
archived LRS downstream workflow instead. It needs the original integrated TSV,
gene/phenotype tables and archived evidence. Neither BAM-free workflow should
be run concurrently with another writer to its output folder.
