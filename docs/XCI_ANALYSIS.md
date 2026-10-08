# X-chromosome inactivation (XCI) analysis

This optional LRS module measures X-chromosome inactivation skew from native-DNA Oxford Nanopore sequencing. It is based on the read-level strategy described by Gocuk et al., Genome Research (2024), "Measuring X-Chromosome inactivation skew for X-linked diseases with adaptive nanopore sequencing" (DOI: 10.1101/gr.279396.124), and on the associated SkewX workflow.

Adaptive sampling is not required for the skew calculation. The published method is applicable to long-read genomic data when there is sufficient chrX coverage. In this repository the analysis uses standard GRCh38 WGS data.

## Biological scope

Run the module only for appropriate XX samples with native-DNA methylation information. The BAM used by the LRS workflow must retain ONT modified-base MM/ML tags. XCI skew is tissue dependent, so blood/saliva/buccal results should not automatically be interpreted as the XCI state of another tissue.

The output is a research measurement of XCI skew. It is not a pathogenicity classification and does not by itself establish whether an X-linked variant is causal.

## Data flow

native-DNA ONT BAM with MM/ML tags
  -> Clair3 SNVs -> WhatsHap phase + haplotag -> read HP/PS labels on chrX
  -> LongPhase SNP/SV phasing -> flip-tolerant phase-consistency QC
  -> CpG-island read-level methylation clustering (NanoMethViz)
       low-methylation epiallele = Xa
       high-methylation epiallele = Xi
       XIST promoter reversed when configured
  -> H1/Xa, H1/Xi, H2/Xa, H2/Xi read counts per phase set
  -> block-wise skew -> folded-binomial maximum-likelihood global skew

The primary estimator is read-level CpG-island methylation linked to read-level haplotypes. Genome-wide/chrX bedMethyl averages are not substituted for this estimator.

## WhatsHap and LongPhase roles

WhatsHap is the primary read-haplotype assignment because the haplotagged BAM provides HP and PS tags for the same reads whose methylation is analyzed.

LongPhase is an independent phase-consistency layer. Haplotype labels are arbitrary between phase blocks and between phasing tools, so agreement is calculated after allowing each overlapping phase-block pair to be either SAME or FLIPPED. A global requirement that H1 from LongPhase equal H1 from WhatsHap would be biologically incorrect.

## Methylation

The XCI branch uses 5mCG for its primary analysis. modkit pileup --phased also generates chrX HP1 and HP2 bedMethyl files from the WhatsHap-haplotagged modBAM. These tracks are used for chromosome-wide visualization/QC.

The read-level estimator uses NanoMethViz clustering over GRCh38 chrX CpG islands. Only islands with exactly two methylation clusters contribute Xa/Xi labels. The default minimum cluster size is 5 reads, matching the published SkewX implementation.

For ordinary CpG islands: lower-methylation cluster -> Xa; higher-methylation cluster -> Xi.
The XIST promoter has the opposite methylation relationship and is reversed when xci_xist_promoter_bed is configured.

## Folded skew

H1_Xa_skew = (H1_Xa + H2_Xi) / (H1_Xa + H1_Xi + H2_Xa + H2_Xi)

Because H1/H2 orientation is arbitrary between phase blocks, the global estimator folds each block around 0.5 and fits the folded-binomial likelihood.

Interpretation:
P = 0.50 -> balanced XCI
P = 0.20 -> approximately 80:20 major:minor XCI
P = 0.10 -> approximately 90:10 major:minor XCI
P approaches 0 -> increasingly extreme skew

The continuous GLOBAL_FOLDED_SKEW_P and XCI_MAJOR_MINOR_RATIO should be reported. Threshold labels in the pipeline are descriptive QC/reporting labels, not clinical decision thresholds.

## Configuration

In snakemake_pipelines/lrs/config_lrs.yaml set xci_samples only for appropriate XX samples and provide a GRCh38 chrX CpG-island BED. xci_xist_promoter_bed is optional but recommended for correct XIST promoter orientation.

## Reference files: CpG islands and XIST promoter

The XCI workflow requires a GRCh38 chrX CpG-island BED. The recommended source is
the UCSC Genome Browser hg38 `cpgIslandExt` table because the workflow follows
the CpG-island methylation strategy used by SkewX. The SkewX publication/workflow
also uses UCSC CpG-island annotation for the X chromosome.

A reproducible server-side download can be made as follows:

```bash
mkdir -p reference/xci

wget \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz \
  -O reference/xci/hg38.cpgIslandExt.txt.gz

gzip -dc reference/xci/hg38.cpgIslandExt.txt.gz \
  | awk 'BEGIN{OFS="\t"} $2=="chrX" {print $2,$3,$4,$5}' \
  | sort -k1,1 -k2,2n \
  > reference/xci/hg38_chrX_cpg_islands.bed
```

The UCSC database dump contains an initial `bin` column, so the BED coordinates
are taken from columns 2-4: chromosome, chromStart and chromEnd. UCSC coordinates
are already 0-based, half-open and therefore compatible with BED format.

Check the result before using it:

```bash
head reference/xci/hg38_chrX_cpg_islands.bed
wc -l reference/xci/hg38_chrX_cpg_islands.bed
```

Then configure:

```yaml
xci_cpg_islands_bed: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/xci/hg38_chrX_cpg_islands.bed"
```

### XIST promoter BED

The XIST promoter BED is optional because there is no single standalone
"official XIST promoter BED" distributed with GRCh38. The purpose of this file
in this workflow is very specific: it marks the XIST promoter CpG island so
that its methylation interpretation is reversed relative to ordinary promoter
CpG islands.

For ordinary X-linked promoter CpG islands:

```text
lower methylation cluster  -> Xa
higher methylation cluster -> Xi
```

For the XIST promoter CpG island:

```text
lower methylation cluster  -> Xi
higher methylation cluster -> Xa
```

XIST is on the minus strand in GRCh38. The MANE/RefSeq locus is
approximately chrX:73820651-73852753, so the transcription start is at the
high-coordinate end. In BED coordinates the one-base MANE TSS can be represented
as:

```text
chrX    73852752    73852753    XIST_TSS
```

The XIST CpG island used for Xa/Xi methylation interpretation does not have to
overlap the exact canonical TSS base. Human XIST studies describe a
differentially methylated CpG island approximately 1.4 kb downstream of the
major P1 promoter, within the 5-prime part of XIST/exon 1. Therefore, testing
only whether a CpG island overlaps one TSS coordinate is too strict and can
legitimately return no result.

Use the UCSC chrX CpG-island BED itself to find CpG islands close to the XIST
5-prime end. XIST is on the minus strand, so transcription proceeds toward
lower genomic coordinates. Search a small window around the GRCh38 XIST
5-prime end and inspect the returned islands:

```bash
awk 'BEGIN{OFS="\t"} \
  $1=="chrX" && $3 > 73847700 && $2 < 73855750 \
  {print $1,$2,$3,$4}' \
  reference/xci/hg38_chrX_cpg_islands.bed
```

The expected XIST-associated island should be close to, but can be downstream
of, the canonical promoter rather than directly overlapping the single TSS
base. Because XIST is on the minus strand, an island about 1.4 kb downstream
will have lower genomic coordinates than the main TSS.

A convenient way to rank nearby CpG islands by distance from the current RefSeq
GRCh38 XIST 5-prime end (approximately 73852714) is:

```bash
awk 'BEGIN{OFS="\t"; tss=73852714}
  $1=="chrX" && $3 > 73847700 && $2 < 73855750 {
    if (tss < $2) d=$2-tss;
    else if (tss > $3) d=tss-$3;
    else d=0;
    print $1,$2,$3,$4,d
  }' reference/xci/hg38_chrX_cpg_islands.bed \
  | sort -k5,5n
```

Do not automatically choose the closest interval only because it is nearest.
Confirm that the selected interval is the CpG island in the 5-prime XIST region
described in the literature and visible in UCSC/GENCODE. Once confirmed, write
that CpG-island interval itself to:

```text
reference/xci/hg38_XIST_promoter_cpg_island.bed
```

For example, after identifying the correct row:

```bash
echo -e "chrX\tSTART\tEND\tXIST_promoter_CpG_island" \
  > reference/xci/hg38_XIST_promoter_cpg_island.bed
```

Replace START and END with the coordinates from the confirmed UCSC CpG-island
row, not with an arbitrary +/-1 kb or +/-2 kb promoter window.

Then configure:

```yaml
xci_xist_promoter_bed: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/xci/hg38_XIST_promoter_cpg_island.bed"
```

If no convincing XIST-associated CpG island is found, keep:

```yaml
xci_xist_promoter_bed: null
```

The workflow still runs without this optional file; only the XIST-specific
methylation reversal is omitted.

## XCI thesis plots

The XCI plotting branch is implemented in `r_plots/plot_xci.R` and is called
by the `xci_plot_analysis` Snakemake rule. Plots are generated only for
samples listed in `xci_samples`.

All XCI figures are written to:

```text
<sample>/plots/15_x_inactivation/
```

Each figure is saved as PDF, SVG and PNG.

### 1. XCI block-skew distribution

```text
<sample>_xci_block_skew_distribution.pdf
```

This is the main summary plot of the XCI estimate. It shows the distribution of
folded block-level skew, weighted by informative read count. The vertical
dotted line at 0.5 represents balanced XCI and the dashed line is the fitted
global `GLOBAL_FOLDED_SKEW_P`.

Interpretation:

```text
P = 0.50 -> approximately 50:50, balanced
P = 0.40 -> approximately 60:40
P = 0.20 -> approximately 80:20
P = 0.10 -> approximately 90:10
P -> 0    -> increasingly strong skew
```

The continuous estimate and `XCI_MAJOR_MINOR_RATIO` should be reported rather
than relying only on the descriptive threshold label.

### 2. XCI along chromosome X

```text
<sample>_xci_chrX_block_skew.pdf
```

This plot places informative phase blocks along chrX. The x-axis is genomic
position in Mb and the y-axis is `H1_Xa_SKEW`.

```text
H1_Xa_SKEW > 0.5 -> H1 is preferentially assigned as Xa
H1_Xa_SKEW < 0.5 -> H2 is preferentially assigned as Xa
H1_Xa_SKEW = 0.5 -> balanced or unresolved orientation
```

Point size reflects the number of informative reads in the block. H1/H2 labels
are local phase labels and may flip between phase blocks.

### 3. Confidence in Xa/Xi orientation

```text
<sample>_xci_orientation_log_odds.pdf
```

This figure plots `LOG10_ODDS_H1_XA_VS_H2_XA` against the number of informative
reads per phase block.

```text
 0  -> equal orientation support
+1  -> at least 10:1 support for H1 as preferential Xa
+2  -> at least 100:1 support for H1 as preferential Xa
-1  -> at least 10:1 support for H2 as preferential Xa
-2  -> at least 100:1 support for H2 as preferential Xa
```

These are orientation likelihoods and must not be interpreted as pathogenicity
probabilities.

### 4. Haplotype-specific chrX methylation

```text
<sample>_xci_haplotype_methylation.pdf
```

This chromosome-wide QC/context plot summarizes CpG 5mC separately for WhatsHap
HP1 and HP2. The default display uses 5-Mb bins from the phased modkit bedMethyl
tracks. It is useful for visualizing broad haplotype-specific methylation
differences but is not the primary XCI skew estimator. The primary estimator
remains read-level CpG-island clustering linked to haplotagged reads.

### 5. WhatsHap versus LongPhase phase concordance

```text
<sample>_xci_phase_concordance.pdf
```

This is a phasing QC plot. It displays the flip-tolerant concordance between
overlapping WhatsHap and LongPhase phase blocks. H1/H2 labels are arbitrary
between independent phasing methods, so phase-block orientation is allowed to
flip before concordance is assessed.

### Recommended thesis use

For the main thesis results, the most informative figures are usually:

```text
xci_block_skew_distribution
xci_chrX_block_skew
xci_haplotype_methylation
```

The orientation-log-odds and WhatsHap/LongPhase concordance plots are useful as
supporting QC figures unless they reveal a biologically important feature.

The XCI result does not by itself identify whether a disease-associated allele
is preferentially active or inactive. For an X-linked candidate SV, the
candidate must first be phased to HP1 or HP2 and then linked to the inferred
Xa/Xi orientation of that haplotype.

## Main outputs

<sample>/xci/
  <sample>_chrX_whatshap_reads.tsv.gz
  <sample>_chrX_whatshap_haplotag_summary.tsv
  <sample>_whatshap_longphase_phase_blocks.tsv
  <sample>_whatshap_longphase_phase_summary.tsv
  <sample>_chrX_clustered_methylation_reads.tsv.gz
  <sample>_chrX_block_skew.raw.tsv.gz
  <sample>_xci_blocks.tsv
  <sample>_xci_summary.tsv
  methylation/<sample>.chrX.hp1.5mC.bedmethyl.gz
  methylation/<sample>.chrX.hp2.5mC.bedmethyl.gz
  methylation/<sample>.chrX.combined.5mC.bedmethyl.gz

Static thesis figures are written under <sample>/plots/15_x_inactivation/ as PDF, SVG and 600-dpi PNG.

## Running

After the upstream LRS update has generated WhatsHap/LongPhase phasing and mosdepth coverage, run:

snakemake --snakefile Snakefile_LRS_postprocess --configfile config_lrs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 xci_analysis

Add -n -p first for a dry run. all_thesis_plots also includes XCI plots for samples listed in xci_samples.

## QC and limitations

The analysis reports chrX coverage, haplotagged-read counts, informative CpG islands/phase blocks, and flip-tolerant WhatsHap/LongPhase concordance.

A sample can legitimately have no estimate when chrX coverage is low, too few reads are haplotagged, CpG islands do not form two interpretable methylation clusters, or too few islands overlap informative phase blocks.

XCI is tissue specific. A skew measured in blood, saliva, or buccal DNA should not automatically be extrapolated to retina, optic nerve, or another tissue.