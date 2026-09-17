# LRS pipeline validation procedure

This checklist validates the repaired long-read workflow **one stage at a time** before running the full cohort. It is intentionally patient-first so failures can be localized without spending compute on every sample.

The repaired data flow is:

```text
BAM
├─ mosdepth
├─ Clair3 → WhatsHap
├─ Sniffles2 discovery ┐
├─ cuteSV genome-wide  ├─ caller evidence filter → Jasmine → MASTER SV VCF
└─ DELLY long-read     ┘                              ├─ AnnotSV
                                                      ├─ VEP
                                                      ├─ LongPhase
                                                      └─ integrated SV/gene table

BAM → Sniffles2 needLR-compatible query → needLR population AF ─┘
BAM → Straglr
BAM → TLDR
BAM → modkit
```

The needLR branch **does not define the master SV universe**. BNDs and SVs >=10 Mb remain available to AnnotSV and the integrated table even though needLR does not evaluate them.

## 1. Work on the repair branch

```bash
git fetch origin
git checkout fix/lrs-pipeline-flow
git pull
cd snakemake_pipelines/lrs
```

## 2. Confirm input files before Snakemake

At minimum, verify:

```bash
test -s /DATA/Reference/hg38.fa
test -s /DATA/Reference/hg38.fa.fai

test -s /DATA/casadei7/tools/SV-PIPELINE-main/reference/human_GRCh38_no_alt_analysis_set.trf.bed
test -s /DATA/casadei7/tools/SV-PIPELINE-main/PANEL_OA/optic_neuropathy_panel.bed
test -s /DATA/casadei7/tools/SV-PIPELINE-main/PANEL_OA/optic_neuropathy_genes.txt

test -d /DATA/casadei7/tools/SV-PIPELINE-main/reference/needlr/release/backend_files
test -d /DATA/casadei7/tools/SV-PIPELINE-main/reference/AnnotSV_annotations/Annotations_Human
```

For every BAM, verify an index exists and that the BAM can be opened:

```bash
samtools quickcheck -v /path/to/sample.bam
samtools idxstats /path/to/sample.bam | head
```

If the BAM index is missing:

```bash
samtools index -@ 16 /path/to/sample.bam
```

## 3. Create/check conda environments

From `snakemake_pipelines/lrs/`:

```bash
snakemake \
  -s Snakefile_LRS_update \
  --use-conda \
  --conda-create-envs-only \
  --cores 1
```

Pay particular attention to:

```text
envs/sniffles2.yaml  → Sniffles2 2.6.2
envs/jasmine.yaml    → Jasmine + Python/pysam + bcftools/tabix
envs/needlr.yaml     → needLR 4.1 + bcftools compatible with needLR 4.1
envs/annotsv.yaml    → AnnotSV + Python + bcftools/htslib/tabix
```

## 4. Full DAG dry-run

```bash
snakemake \
  -s Snakefile_LRS_update \
  --use-conda \
  --cores 1 \
  --dry-run \
  --printshellcmds \
  --reason
```

Do not proceed to the full cohort if Snakemake reports ambiguous rules, missing inputs, duplicate outputs, wildcard errors or environment-resolution errors.

## 5. Test one patient through caller parsing/filtering

Using sample `14996` from the current config:

```bash
ROOT=/DATA/casadei7/tools/SV-PIPELINE-main/outs_test_CNTRL
SAMPLE=14996
```

Run only the three filtered caller VCFs:

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --printshellcmds \
  "$ROOT/$SAMPLE/sv/sniffles2/${SAMPLE}_sniffles.filtered.vcf" \
  "$ROOT/$SAMPLE/sv/cutesv/${SAMPLE}_cutesv.filtered.vcf" \
  "$ROOT/$SAMPLE/sv/delly/${SAMPLE}_delly.filtered.vcf"
```

Check that every caller has unique, non-missing IDs:

```bash
for VCF in \
  "$ROOT/$SAMPLE/sv/sniffles2/${SAMPLE}_sniffles.vcf" \
  "$ROOT/$SAMPLE/sv/cutesv/${SAMPLE}_cutesv.vcf" \
  "$ROOT/$SAMPLE/sv/delly/${SAMPLE}_delly.vcf"
do
  echo "=== $VCF ==="
  bcftools query -f '%ID\n' "$VCF" | awk '$1=="." || $1=="" {n++} END{print "missing IDs:",n+0}'
  bcftools query -f '%ID\n' "$VCF" | sort | uniq -d | head
 done
```

Inspect PASS/FAIL reasons:

```bash
column -t -s $'\t' "$ROOT/$SAMPLE/parsed/${SAMPLE}.sniffles.evidence.filtered.tsv" | head
column -t -s $'\t' "$ROOT/$SAMPLE/parsed/${SAMPLE}.cutesv.evidence.filtered.tsv" | head
column -t -s $'\t' "$ROOT/$SAMPLE/parsed/${SAMPLE}.delly.evidence.filtered.tsv" | head
```

The DELLY table should now obtain alternate support from DELLY evidence (e.g. FORMAT/DV when INFO/SU is absent), not from reference depth.

## 6. Test Jasmine independently

Run the master merged VCF:

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --printshellcmds \
  "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_merged_SV.vcf.gz"
```

Validate compression/index/header/order:

```bash
bcftools view -h "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_merged_SV.vcf.gz" | tail
bcftools index -n "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_merged_SV.vcf.gz"
bcftools stats "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_merged_SV.vcf.gz" > "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_merged_SV.bcftools.stats.txt"
```

Then generate caller concordance:

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 4 --printshellcmds \
  "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_caller_support_summary.tsv"
```

Inspect `SUPP`, `SUPP_VEC`, decoded `CALLERS`, and `CALLER_COUNT`:

```bash
column -t -s $'\t' "$ROOT/$SAMPLE/sv/merged/${SAMPLE}_caller_support_summary.tsv" | head -30
```

The configured `jasmine_caller_order` must match the file-list order exactly:

```text
Sniffles2
cuteSV
delly
```

The `high_confidence` VCF is a **companion review set only**. Do not use it as the only clinical/discovery callset.

## 7. Test the dedicated needLR branch

First produce the dedicated Sniffles2 query VCF:

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --printshellcmds \
  "$ROOT/$SAMPLE/sv/needlr/${SAMPLE}_needLR_query.sniffles.vcf"
```

Then run needLR:

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --printshellcmds \
  "$ROOT/$SAMPLE/sv/needlr/${SAMPLE}_needLR_RESULTS.tsv"
```

Inspect the TSV header:

```bash
head -2 "$ROOT/$SAMPLE/sv/needlr/${SAMPLE}_needLR_RESULTS.tsv" | column -t -s $'\t'
```

Expected current needLR-style fields include coordinates/type/length, query ID and an all-population AF field such as `Allele_Freq_ALL`.

Important interpretation:

- BNDs are not evaluated by needLR.
- SVs >=10 Mb are not evaluated by needLR.
- A missing needLR match must **not** automatically be interpreted as AF=0.
- Non-zero chrX/chrY AFs receive a caution in the integrated table because of needLR's denominator limitation.

## 8. Test AnnotSV directly from the Jasmine master VCF

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 16 --printshellcmds \
  "$ROOT/$SAMPLE/sv/annotsv/${SAMPLE}_merged_SV.annotsv.tsv"
```

This rule must depend on:

```text
.../sv/merged/<sample>_merged_SV.vcf.gz
```

and **not** on the needLR result VCF.

Inspect the header and gene column:

```bash
head -2 "$ROOT/$SAMPLE/sv/annotsv/${SAMPLE}_merged_SV.annotsv.tsv" | cut -c1-500
```

## 9. Test the integrated table

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 16 --printshellcmds \
  "$ROOT/$SAMPLE/gene_discovery/${SAMPLE}_integrated_SV_gene_analysis.tsv"
```

Inspect key columns:

```bash
head -1 "$ROOT/$SAMPLE/gene_discovery/${SAMPLE}_integrated_SV_gene_analysis.tsv" | tr '\t' '\n'
```

Key provenance fields should include:

```text
SV_ID
CHROM
START
END
SVTYPE
SVLEN
CALLERS
CALLER_COUNT
CALLER_SUPPORT_CLASS
SUPP_VEC
EXACT_ID_READ_SUPPORT
GENES
NEEDLR_STATUS
NEEDLR_AF
NEEDLR_QUERY_ID
NEEDLR_MATCH_METHOD
NEEDLR_BP_DISTANCE
NEEDLR_REL_SIZE_DIFFERENCE
NEEDLR_AF_CAUTION
PANEL_STATUS
PHENOTYPE_SCORE
CANDIDATE_CLASS
```

Check needLR matching status counts:

```bash
awk -F'\t' '
NR==1 {for(i=1;i<=NF;i++) if($i=="NEEDLR_STATUS") c=i; next}
{n[$c]++}
END {for(k in n) print k,n[k]}
' "$ROOT/$SAMPLE/gene_discovery/${SAMPLE}_integrated_SV_gene_analysis.tsv" | sort
```

`NO_MATCH` is different from `AF=0` and must remain distinct.

## 10. Test the other biological branches

```bash
snakemake -s Snakefile_LRS_update --use-conda --cores 32 --printshellcmds \
  "$ROOT/$SAMPLE/sv/straglr/${SAMPLE}_straglr.annotated.tsv" \
  "$ROOT/$SAMPLE/phasing_longphase/${SAMPLE}.longphase.vcf.gz" \
  "$ROOT/$SAMPLE/mei/tldr/${SAMPLE}.tldr.table.txt" \
  "$ROOT/$SAMPLE/methylation/${SAMPLE}.modkit.tsv"
```

The Straglr annotation script now parses the variable-width gene BED after the four-column locus BED instead of using fixed incorrect intersection offsets.

## 11. Only after one patient passes, run the full configured set

```bash
snakemake \
  -s Snakefile_LRS_update \
  --use-conda \
  --cores 32 \
  --printshellcmds \
  --rerun-incomplete \
  --keep-going
```

For debugging, omit `--keep-going` so the first failing rule stops the workflow.

## 12. Preserve the original results

The repaired workflow writes the same major canonical outputs but never requires deleting the original `Snakefile_LRS`. Keep previous result directories separate while validating this branch. If possible, point `path:` in `config_lrs.yaml` to a fresh test directory for the first validation run.
