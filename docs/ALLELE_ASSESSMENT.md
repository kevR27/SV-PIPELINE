# Allele assessment: configuration, interpretation and migration

For output-only studies and the distinction between gene ranking, AnnotSV CNV
classification and the additional allele score, see [Score interpretation](SCORE_INTERPRETATION.md).

The original question, ONT/SRS comparison, genome-wide SV discovery, human-only
phenotype evidence and optic-neuropathy panel remain intact. This addition
assesses an SV–gene hypothesis rather than replacing the master callset or
claiming clinical pathogenicity. It retains all rows, including intergenic events,
BNDs, large SVs and outside-panel genes. VEP retains `--flag_pick`; AnnotSV remains
the genome-wide gene mapper. No additional caller or non-human data is required.

## Outputs

In `<output-root>/<sample>/gene_discovery/`:

- `<sample>_integrated_SV_gene_analysis.tsv`: original integration, with corrected
  joins, source fields and structured caller/transcript evidence.
- `<sample>_allele_assessment.tsv`: same row order/count plus the selected hypothesis,
  six evidence domains, score, dense research rank, missing domains and review flags.
- `<sample>_allele_disease_hypotheses.tsv`: one row per SV–gene–curated disease model,
  ordered by score. No mixing of different diseases' inheritance/mechanism evidence.
  Unknown-disease rows remain if no model exists. The highest-scoring hypothesis is
  selected for the main table; ties use disease ID, MOI and mechanism deterministically.
- `<sample>_allele_assessment.manifest.json`: model version, thresholds, input paths
  and SHA256 hashes. Preserve with results and record the repository commit.

LRS postprocessing and SRS complementary-evidence integration consume the new
assessment. Existing gene-discovery scores remain separate. A plot of the old
`ranked_candidates.tsv` is still a gene-ranking plot, not the new allele ranking.
SV–gene rows, disease-hypothesis rows, and paired BND records are not unique molecular
event counts. Event-level adjudication/deduplication is still required for thesis totals.

## Research score v1

Each domain contributes 0–2 points; the sum has a fixed maximum of 12. It is not
normalized by data availability. Unknown evidence contributes no points and is
explicitly listed, rather than treated as negative or evidence of rarity. Evidence
categories can be correlated: this is an uncalibrated review aid, not a probability,
clinical classification, or automatic exclusion rule. Evaluate thresholds and ranking
on known positives/held-out cases before relying on rank.

| Domain | Assessment and limits |
|---|---|
| Population | An external source-backed exact/sequence-validated allele AF ≤ threshold gives 2; zero observed AF gives 1, with an explicit absence-not-pathogenicity label. AF above threshold gives 0 and a flag. needLR coordinate-only AF remains visible as provisional and gives 0 until allele identity is reviewed. Missing/unmatched/ambiguous is unknown. A benign database overlap AF is not transferred as the patient's exact allele AF. |
| Technical | At least the configured maximum per-caller read support gives 1 for a single caller or 2 with ≥2 callers. Low support, ambiguity and QC flags receive 0 and a review status. Read counts across callers are not summed. Agreement on the same reads is not independent validation. These thresholds do not change caller/filter thresholds. |
| Disruption | AnnotSV gene-specific transcript rows can identify CDS loss, transcript copy gain or possible exonic breakpoint disruption, giving 1 as predictions. A source-backed reviewed effect gives 1 if predicted or 2 if validated; no disruption gives 0. Gene overlap alone does not prove LOF. Partial duplications, inversions, insertions, regulatory effects and GOF/dominant-negative mechanisms often need manual molecular review. |
| Inheritance | QC-passing exact caller-linked or externally supplied genotypes are evaluated against one curated disease MOI. AD genotype compatibility gives 1; callable reference parents produce a candidate de novo flag with 2, not confirmed de novo status. AR homozygous-alt compatibility gives 2; a single heterozygous variant requires another allele. Two distinct nonoverlapping same-gene SVs with compatible effects and opposite QC-passing parental transmissions are flagged as inferred trans candidates (1), not confirmed causality. |
| Patient phenotype | Deduplicated patient positive HPO terms are compared with disease-specific HPO terms, otherwise human gene associations. ≥50% exact matching gives 2; any smaller positive match gives 1. A matching explicitly absent term gives a review flag and 0. No exact match is not proof of incompatibility. This v1 implementation is exact-term matching, not ontology-aware semantic similarity. Generic optic anchors are never substituted for patient phenotypes. |
| Disease mechanism | Source-backed disease models with compatible effect and Moderate/Strong/Definitive validity give 2; Limited validity or a possible effect gives 1. Disputed/refuted associations or mechanism mismatch are flagged and give 0. ClinGen dosage evidence without a curated disease assignment can give 1, explicitly labeled gene-level only. |

Genotype QC defaults: GQ ≥20, depth ≥5. Missing QC is unknown. Merged Jasmine/SURVIVOR
genotypes are never used as a patient genotype. Coordinate-only caller links do not
transfer genotype. Conflicting high-quality caller genotypes remain conflict.
Parent absence from a VCF is not a reference genotype. Parent genotypes require exact
master SV IDs, callable reference/alternate genotyping, source and QC.

Sex-linked/PAR/sex-chromosome dosage and mitochondrial inheritance remain review
categories; ordinary diploid rules are not applied automatically. Family identity,
parentage, mosaicism and independent validation must still be assessed. SNV–SV
compound heterozygosity and direct LongPhase phase-block inference are not automated
by this stage. LongPhase/WhatsHap outputs remain complementary context downstream.

## Configuration

The same keys are supported in `config_lrs.yaml` and `config_srs.yaml`:

```yaml
allele_patient_context: /absolute/path/patient_context.tsv
allele_disease_models: /absolute/path/disease_models.tsv
allele_evidence: null
allele_family_genotypes: null
allele_rare_af: 0.01
allele_min_support: 5
allele_min_gq: 20
allele_min_dp: 5
```

All four input files are optional. `null` runs the assessment with explicit unknowns.
An optional path can instead be a sample-ID mapping, allowing files to be supplied
only for selected samples. Samples not in that mapping have that evidence unavailable.
A configured patient-context file must contain the requested sample. Missing configured
files are errors, not silently converted to missing evidence. Use the header-only TSV
templates under `reference/allele_assessment/`; fill with your actual reviewed data.
Do not put the template paths in config until they contain appropriate patient rows.

### Patient context

Columns: `sample_id`, `hpo_present`, `hpo_absent`, `sex`, `mother_id`, `father_id`, `source`.
Use exact sample identifiers and semicolon-separated HP:NNNNNNN terms. `.` means
unknown/not supplied. Absence must be explicitly assessed; do not mark every unrecorded
phenotype absent. Provide record/source provenance. The same term cannot be both present
and absent. All associations remain human-only.

### Disease models

Columns: `gene`, `disease_id`, `moi`, `mechanism`, `validity`, `hpo_terms`, `source`.
One reviewed gene–disease–MOI–mechanism hypothesis per row. Use an OMIM/MONDO or other
stable disease ID and a versioned/source citation. Allowed values:

- `moi`: AD, AR, XLD, XLR, MT, UNKNOWN.
- `mechanism`: LOF, HAPLOINSUFFICIENCY, TRIPLOSENSITIVITY, GOF, DOMINANT_NEGATIVE, OTHER, UNKNOWN.
- `validity`: DEFINITIVE, STRONG, MODERATE, LIMITED, DISPUTED, REFUTED, UNKNOWN.

HPO terms are optional. A missing disease model is not a negative gene–disease result.
AnnotSV's aggregate GenCC/OMIM strings are retained for review but are not automatically
zipped together into disease-specific assertions: multi-disease lists do not reliably
establish the linkage of disease, MOI and molecular mechanism. Curate/import these rows
from source records that preserve those relationships. This is deliberately a reviewed
input rather than a fabricated automatic mapping.

### Allele evidence

Columns: `sample_id`, `SV_ID`, `gene`, `disease_id`, `effect`, `effect_status`, `gt`,
`gq`, `dp`, `population_af`, `population_match`, `population_source`, `source`.
An exact disease row overrides the same allele/gene row with disease_id `.`; rows are
not merged field-by-field. Duplicate keys are errors. Gene symbols are normalized to
uppercase; SV IDs and sample IDs are exact.

`effect`: LOF, COPY_GAIN, GOF, DOMINANT_NEGATIVE, NO_DISRUPTION, OTHER, UNKNOWN (or `.`).
A specified effect requires `effect_status` PREDICTED or VALIDATED and a real source.
A validated effect means evidence for that molecular effect, not merely PCR confirmation
that a structural event exists. Optional genotype uses VCF GT notation and GQ/depth.
AF is a number from 0 to 1. An AF contribution requires `population_match` EXACT_ALLELE
or SEQUENCE_VALIDATED and a population source. Record the relevant dataset release,
ancestry/denominator/callability limitations in the source/review record; the script does
not independently verify external assertions or independently query databases.

### Family genotypes

Columns: `sample_id`, `SV_ID`, `gt`, `gq`, `dp`, `source`.
Supply child/parent genotype records produced by an appropriate joint/targeted genotyping
or reviewed workflow against the same master allele identity. Absence of a record is
unknown. These records take precedence over caller genotype summaries, including when
they fail QC. Do not fill missing genotypes with 0/0.

## Correctness fixes and migration

- Exact panel membership prevents `NONPANEL_GENE` matching `PANEL_GENE`.
- AnnotSV full records supply event-level fields; only exact single-gene records supply
  gene-level fields. All gene transcript records are retained as JSON. Conflicting scalar
  transcript values remain missing; the disruption assessment inspects the full records.
- BND parsing reads ALT partners/orientation; coordinate fallback requires both endpoints
  and orientation. Unknown orientation is not assumed. Ambiguous caller matches are not
  arbitrarily assigned. Source caller JSON includes read support, genotype/QC and IDs. Jasmine prefixed
  IDs are resolved using the explicit caller file-list order; do not reorder the
  `--caller-tsv` arguments relative to Jasmine input files.
- needLR matching remains coordinate/size based. Multiple candidates become ambiguous;
  a unique match is still provisional, particularly for insertion sequence identity.
  The same AF parser is used for integration and plots, with original field provenance.
  Invalid frequencies never become zero. Match failures never imply population absence.
- Duplicate gene–HPO associations no longer inflate the original gene score.
- Gene summaries report database YES, NOT_REPORTED, UNKNOWN and NOT_APPLICABLE counts,
  plus unrecognized legacy values. `NOT_REPORTED` is not an assessed negative search.
- LongPhase declares/indexes both SNP and SV VCFs. SV matching/plots use `_SV.vcf.gz`
  (or the legacy uncompressed `_SV.vcf`), never the SNP output as a fallback.
- modkit bedMethyl units are explicitly `percent` by default, in integration and plots.
  Use `methylation_units: fraction` only for a genuinely converted fractional input.
  BND methylation/nearby-phasing context uses the correct chromosome for each endpoint.
  Repeat/MEI contextual overlaps for BNDs are explicitly first-breakpoint context,
  not a chromosome-spanning interval or whole-rearrangement confirmation.
- needLR executes in a fresh retained run directory whenever its rule is rerun. It no
  longer trusts an arbitrary old native result. Backend files are declared inputs; set
  `needlr_backend_version` to the installed release. The native run/query checksum is
  recorded. Old run directories remain and may consume disk space; clean them manually
  only after reviewing what you need to retain. Backend checksums are not generated;
  retain a versioned backend manifest alongside the configured release.
- Optional LongPhase, methylation, index and WhatsHap files are tracked as workflow inputs.
  The postprocessing workflow uses only files already present at DAG construction time;
  complete the main pipeline before running it.
- `Snakefile_LRS` is now a compatibility include of `Snakefile_LRS_update`, avoiding two
  divergent implementations. The canonical workflow's existing configuration is required.

After updating, rerun the integrated and assessment targets to regenerate affected
results, then postprocessing. If resuming old outputs, force the changed normalization,
integration and ranking rules where necessary; do not assume an old table has the new
schema merely because its filename is unchanged. This may cause needLR to rerun once
because the rule/source changed. Main calling thresholds and discovery scope are unchanged.

Example, from `snakemake_pipelines/lrs/`, using your configured sample/output paths:

```bash
snakemake -s Snakefile_LRS_update --configfile config_lrs.yaml \
  --use-conda --cores 32 /absolute/output-root/SAMPLE/gene_discovery/SAMPLE_allele_assessment.tsv
snakemake -s Snakefile_LRS_postprocess --configfile config_lrs.yaml \
  --use-conda --cores 32 all_thesis_plots
```

The assessment can also run directly using `python scripts/assess_sv_alleles.py --help`.
Use the same reviewed configuration and optional inputs for both platforms. Keep
configuration files with patient paths/data out of public commits.

## Verification and remaining research work

Run `python -m unittest discover -s tests -v` with numpy/pandas/matplotlib installed.
Fixtures cover missingness, ambiguity, multi-gene scope, BND identity, duplicate HPOs,
units, genotype conflicts, inheritance, phenotype and end-to-end output preservation.
The CI also parses both main workflows and postprocessing with Snakemake 8.30.0.
These checks do not validate patient biology, SV calling accuracy, the ranking model,
or a full server run. Benchmark known positives, use matched LRS/SRS cases, and validate
key findings with the study team before drawing diagnostic conclusions.

Primary references: [ACMG/ClinGen CNV standards](https://doi.org/10.1038/s41436-019-0686-8),
[ClinGen dosage framework](https://www.clinicalgenome.org/curation-activities/dosage-sensitivity/),
[AnnotSV documentation](https://lbgi.fr/AnnotSV/),
[LongPhase output source](https://github.com/twolinin/longphase/blob/v1.7.3/ParsingBam.cpp),
[modkit bedMethyl specification](https://github.com/nanoporetech/modkit/blob/master/book/src/intro_pileup.md).
