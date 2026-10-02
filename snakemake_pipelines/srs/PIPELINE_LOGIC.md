# SRS Illumina WGS structural-variant workflow

This document describes the **SRS-only** workflow in `Snakefile_SRS`.
The LRS workflows are independent and are not imported or modified by this
implementation.

## Biological objective

The SRS workflow is designed for rare-disease WGS in hereditary optic
neuropathy. It keeps the same research logic used for the LRS analysis:

1. preserve a genome-wide discovery universe;
2. distinguish known optic-neuropathy panel genes from non-panel discovery;
3. retain large and complex events rather than silently filtering them;
4. separate technical support from biological interpretation;
5. combine independent evidence modalities when they are genuinely
   complementary;
6. retain mitochondrial and nuclear-mitochondrial analyses;
7. preserve unmatched findings from complementary callers instead of forcing
   every tool into one merged VCF.

It is **not** an attempt to recreate GATK-SV and does not require GATK-SV.

## Source-of-truth hierarchy

### 1. Nuclear breakpoint master: Manta + DELLY -> Jasmine

Manta and DELLY are the primary short-read breakpoint callers. Their PASS
calls are normalized into a common evidence table and restricted to canonical
chromosomes (chr1-22, chrX, chrY, chrM).

Pre-merge support is intentionally permissive: a PASS call with at least one
explicit supporting read signal can remain in discovery. Low-support calls are
not automatically interpreted as reliable; downstream evidence review remains
necessary.

Jasmine merges the two callers using size-aware breakpoint distance:

- `max_dist = 5000`
- `max_dist_linear = 0.5`
- `min_dist = 50`
- `--mutual_distance`
- `--normalize_type`
- `--allow_intrasample`

The complete Jasmine union is the row-defining master VCF.

The companion `*.caller_concordant.vcf.gz` contains events supported by both
Manta and DELLY. It is intentionally **not** called "high confidence":
agreement between two callers that use overlapping PE/SR evidence is not an
independent validation assay.

References:

- Manta: https://github.com/Illumina/manta
- DELLY: https://github.com/dellytools/delly
- Jasmine: https://github.com/mkirsche/Jasmine

### 2. Functional annotation: AnnotSV + VEP

AnnotSV remains the primary genome-wide SV/gene annotation layer.

VEP is supplementary and is run without collapsing the output to a single
picked transcript. Transcript consequences, exon/intron context, canonical
status and VEP PICK flags are retained and summarized per SV-gene pair by
`summarize_vep_sv_consequences.py`.

VEP is therefore not used as the only source of affected genes.

### 3. Independent CNV evidence: CNVpytor

CNVpytor adds read-depth evidence that Manta/DELLY cannot replace. The default
workflow uses:

- 1 kb read-depth calling;
- genome-wide DeepVariant PASS biallelic SNPs as BAF input;
- 10 kb BAF evaluation/genotyping.

CNVpytor raw calls are retained. Q0, pN and distance-to-reference-gap are
review flags, not hard pathogenicity or validity filters.

CNVpytor DEL/DUP calls are matched back to the Jasmine master only when their
geometry is compatible. CNVpytor calls with no master match are written to:

`gene_discovery/<sample>_cnvpytor_only_candidates.tsv`

and are annotated separately with AnnotSV.

This prevents a real read-depth CNV from disappearing merely because a
breakpoint caller did not find it.

Reference:

- CNVpytor: Suvakov et al., GigaScience 2021, DOI 10.1093/gigascience/giab074
- Project: https://github.com/abyzovlab/CNVpytor

### 4. Independent breakpoint-assembly evidence: GRIDSS

GRIDSS is deliberately **not** inserted as a third equivalent member of the
Manta/DELLY Jasmine merge.

GRIDSS natively represents rearrangements as breakends and uses local assembly.
Forcing every GRIDSS BND into DEL/DUP/INV equivalence would add an
interpretation step that is not justified by the caller output.

The workflow therefore:

1. runs GRIDSS independently;
2. collapses reciprocal GRIDSS breakend records into descriptive event rows;
3. attaches a GRIDSS event to a master event only when breakpoint geometry is
   compatible;
4. writes unmatched GRIDSS events to
   `gene_discovery/<sample>_gridss_only_events.tsv`.

GRIDSS AS/RAS/SR/RP/VF and FILTER/QUAL evidence are preserved. No new
home-made GRIDSS confidence class is generated.

Reference:

- GRIDSS: https://github.com/PapenfussLab/gridss

### 5. Small variants and BAF: genome-wide DeepVariant

DeepVariant is run genome-wide with the WGS model. The optic-neuropathy panel
VCF is derived **after** genome-wide calling.

This avoids the previous asymmetry in which SV discovery was genome-wide but
SNV/indel calling was panel-restricted.

The WGS DeepVariant callset is used for:

- the complete SNV/indel callset;
- panel extraction;
- CNVpytor BAF input;
- short-read WhatsHap phasing.

Short-read phasing is supportive local evidence. It must not be interpreted as
equivalent to long-read haplotype resolution.

### 6. mtDNA: mity

mity is used as a dedicated mitochondrial WGS SNV/indel and heteroplasmy
branch. It is explicitly configured for `hg38`; GRCh38 chrM is rCRS
(NC_012920.1, 16,569 bp).

The pipeline preserves:

- mity FILTER;
- mity VAF rather than redefining heteroplasmy;
- tier and q;
- mity position/strand/mapping/base-quality filter fields;
- the native mity Excel report.

The mtDNA callset stays separate from the nuclear structural-variant table.
It should be interpreted as a mitochondrial result set, not forced into
nuclear SV semantics.

Reference:

- mity: https://github.com/KCCG/mity

### 7. Specialized branches retained

- **MELTv2**: mobile-element insertions. Disabled until manually installed
  resources are configured.
- **ExpansionHunter**: targeted genotyping of loci in the configured repeat
  catalog. It is not represented as de novo genome-wide repeat discovery.
- **MitoCarta 3.0**: optional context for nuclear genes affecting mitochondrial
  biology. It is disabled until local MitoCarta files are configured.

## Final integrated SRS evidence

The main nuclear interpretation output is:

`gene_discovery/<sample>_integrated_SRS_evidence.tsv`

The master Jasmine SV remains row-defining. Added SRS evidence includes:

- Manta/DELLY caller provenance;
- caller read support and flags;
- AnnotSV database/gene evidence;
- VEP transcript context;
- panel/non-panel status;
- phenotype/gene-discovery context;
- CNVpytor read-depth evidence;
- CNVpytor BAF context;
- GRIDSS local-assembly/breakpoint support;
- MELT/ExpansionHunter contextual matches where applicable.

`SRS_EVIDENCE_MODALITIES` is descriptive. Its count is **not** a
pathogenicity score and is not a validated clinical classifier.

## Internal cohort recurrence

The workflow also creates:

- `cohort/srs_sv_recurrence.tsv`
- `cohort/srs_sv_recurrence.members.tsv`

This is a conservative within-study recurrence analysis. It is useful for
recognizing recurrent candidates and possible technical artifacts.

It is **not** population allele frequency and must not be substituted for
gnomAD-SV/DGV/other population evidence.

## WGS QC

The workflow records:

- mosdepth mean coverage;
- median depth estimated from the cumulative global distribution;
- fraction and percentage of bases >=10x, >=20x and >=30x;
- total/mapped/properly-paired reads;
- duplicates;
- mean/SD insert size;
- samtools error rate;
- idxstats.

The summary is descriptive and does not invent universal PASS/FAIL clinical
thresholds.

Somalier/VerifyBamID-style sample identity, ancestry, sex and contamination QC
are not fabricated here because they require the appropriate validated
reference/site resources. Existing Somalier cohort QC can remain a separate
pre-analysis QC layer until those resources are explicitly configured for this
workflow.

## Important non-claims

The workflow does not infer:

- pathogenicity from the number of callers;
- rarity from absence in the local cohort;
- inheritance when parental/trio data are unavailable;
- cis/trans relationships that short-read phase blocks do not resolve;
- direct gene disruption from a large inversion merely because a gene lies
  inside the inverted interval;
- exact-allele population equivalence from a loose database overlap.

These distinctions should be retained in thesis Methods and Discussion.

## Required server resources

Before a production run, verify:

1. coordinate-sorted BAM and matching BAI for every sample;
2. GRCh38/hg38 FASTA and FAI;
3. BWA indexes next to the FASTA for GRIDSS:
   `.amb .ann .bwt .pac .sa`;
4. AnnotSV annotations;
5. offline VEP GRCh38 cache;
6. Monarch nodes/edges used by the existing phenotype layer;
7. DELLY GRCh38 exclusion BED;
8. ExpansionHunter catalog;
9. optional MELTv2 resources if MELT is enabled;
10. optional MitoCarta workbook/GMX files if MitoCarta is enabled.

The Bioconda CNVpytor recipe includes its reference-resource data, so the
workflow does not invoke `cnvpytor -download` during a run.

## Recommended first server test

Do not begin with the full cohort. Configure one representative Illumina WGS
sample and run:

```bash
snakemake \
  -s Snakefile_SRS \
  --configfile config_srs.yaml \
  --use-conda \
  --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32 \
  all
```

For the first sample, inspect at minimum:

- `qc/<sample>.srs_wgs_qc.tsv`
- `sv/merged/<sample>_caller_support_summary.tsv`
- `cnv/cnvpytor/<sample>.cnvpytor.calls.tsv`
- `gene_discovery/<sample>_cnvpytor_only_candidates.tsv`
- `sv/gridss/<sample>.gridss.events.tsv`
- `gene_discovery/<sample>_gridss_only_events.tsv`
- `mtdna/mity/<sample>.mity.tsv`
- `mtdna/mity/<sample>.mity.report.xlsx`
- `gene_discovery/<sample>_integrated_SRS_evidence.tsv`

Only after the output distributions and candidate examples are biologically
plausible should the workflow be expanded to the complete cohort.
