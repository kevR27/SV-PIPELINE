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

The XIST-specific interval is optional for running the workflow, but it is
biologically useful because XIST has the opposite methylation relationship from
ordinary X-linked promoter CpG islands.

For ordinary promoter CpG islands:

```text
lower methylation cluster  -> Xa
higher methylation cluster -> Xi
```

For the XIST P2 differentially methylated CpG island:

```text
lower methylation cluster  -> Xi
higher methylation cluster -> Xa
```

This XIST interval is not present in the standard UCSC `cpgIslandExt` chrX
subset because the published XIST P2 region is an intermediate-density CpG
island defined with more relaxed CpG-island criteria.

Chapman et al. (2014) show the human XIST P2 CpG-island boundaries in hg19 as
approximately chrX:73,070,847-73,071,444 in genomic position coordinates.
Because BED uses a 0-based start, first represent the published hg19 interval as:

```text
chrX    73070846    73071444    XIST_P2_intermediate_CpG_island
```

Do not guess the GRCh38 coordinates manually for production use. Convert this
published hg19 interval with UCSC liftOver:

```bash
mkdir -p reference/xci/liftover
cd reference/xci/liftover

wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

wget \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

printf "chrX\t73070846\t73071444\tXIST_P2_intermediate_CpG_island\n" \
  > XIST_P2.hg19.bed

./liftOver \
  XIST_P2.hg19.bed \
  hg19ToHg38.over.chain.gz \
  XIST_P2.hg38.bed \
  XIST_P2.unmapped.bed

cat XIST_P2.hg38.bed
cat XIST_P2.unmapped.bed
```

The GRCh37 and GRCh38 RefSeq XIST loci differ by an overall +780,165 bp offset
at both annotated ends, so an unchanged local mapping would be expected near
GRCh38 chrX:73,851,012-73,851,609 in 1-based genomic coordinates
(approximately BED chrX:73851011-73851609). The liftOver result is authoritative
and should be used instead of the expected coordinates.

Copy the successful lifted row to:

```text
reference/xci/hg38_XIST_P2_CpG_island.bed
```

The standard UCSC chrX CpG-island input does not contain this intermediate-density
P2 island. Therefore the lifted XIST interval must also be appended to the CpG
clustering input; otherwise `xci_xist_promoter_bed` has nothing to reverse:

```bash
cat \
  reference/xci/hg38_chrX_cpg_islands.bed \
  reference/xci/hg38_XIST_P2_CpG_island.bed \
  | sort -k1,1 -k2,2n -k3,3n \
  | uniq \
  > reference/xci/hg38_chrX_cpg_islands.with_XIST_P2.bed
```

Configure both files:

```yaml
xci_cpg_islands_bed: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/xci/hg38_chrX_cpg_islands.with_XIST_P2.bed"
xci_xist_promoter_bed: "/DATA/casadei7/tools/SV-PIPELINE-main_v3/reference/xci/hg38_XIST_P2_CpG_island.bed"
```

If `xci_xist_promoter_bed` is left null, the global XCI analysis still runs.
The missing file affects only the XIST-specific reversed methylation assignment;
it does not disable haplotagging, methylation clustering at the other chrX CpG
islands, block-wise skew calculation, the folded-binomial global estimate, or
the other XCI QC outputs.

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