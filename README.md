# SV-PIPELINE: Whole-genome variant analysis

This repository contains two Snakemake workflows: **long-read sequencing (LRS)** from Oxford Nanopore and **short-read sequencing (SRS)** from Illumina. They analyze whole-genome data, identify structural variants (SVs), and organize evidence for hereditary optic neuropathy (HON) research. Nuclear genes involved in mitochondrial biology are an important research focus; mitochondrial-DNA findings are reviewed separately.

Variant discovery is genome-wide. The HON gene panel is used after discovery to label and prioritize genes, not to limit the genomic regions analyzed. Results are research candidates and still require review and, where appropriate, independent confirmation.

This page explains the LRS workflow and the shared code. For Illumina-specific steps, see [the SRS README](snakemake_pipelines/srs_wgs_pipeline/README.md). For detailed figure descriptions, see [the plotting guide](plots/README.md).

## 1. Where the files are

| Location | Purpose |
| --- | --- |
| **snakemake_pipelines/lrs/Snakefile_LRS_update** | Main ONT analysis: QC, variant calling, merging and annotation. |
| **snakemake_pipelines/lrs/Snakefile_LRS_postprocess** | Combines evidence, ranks candidates and prepares downstream results. |
| **snakemake_pipelines/lrs/Snakefile_LRS** | Short compatibility entry point that uses the active LRS workflow. |
| **snakemake_pipelines/lrs/config_lrs.yaml** | Sample BAMs, output folder, reference files, thresholds and optional analyses. |
| **snakemake_pipelines/lrs/thesis_extensions.smk** | Additional depth, HPO, SNV–SV and final-ranking rules. |
| **snakemake_pipelines/lrs/xci_analysis.smk** | Chromosome-X phasing, methylation and XCI analysis. |
| **snakemake_pipelines/srs_wgs_pipeline/** | Separate Illumina workflow and configuration. |
| **snakemake_pipelines/common/** | Rules shared by both workflows: allele assessment, AnnotSV auditing, Exomiser and optional AlphaGenome. |
| **scripts/** | Python and R scripts that read, filter, annotate and summarize the results. |
| **plots/** and **r_plots/** | Scripts for figures, candidate review and visual summaries. |
| **envs/** | Conda environment definitions for individual tools. |
| **PANEL_OA/** and **reference/** | HON panel files and reference/annotation resources. |
| **docs/** | Additional detailed guides when manual setup or interpretation needs more explanation. |

In a Snakefile, **input** is what a step reads, **output** is what it writes, **threads** is its CPU allocation, **conda** selects its software environment, and **shell** gives the executed command. Snakemake uses file dependencies to determine the order automatically.

## 2. Inputs and configuration

The LRS workflow starts with **coordinate-sorted, indexed ONT BAM files** aligned to GRCh38. Edit **config_lrs.yaml** to set the sample-to-BAM mapping, result folder (**path**), reference FASTA (**ref**), available CPU threads, and annotation/reference paths. Quote numeric sample IDs, such as "21930". BAM indexes and reference files must match the alignments. For methylation and XCI analyses, the BAM also needs the appropriate modified-base information.

External reference data, such as AnnotSV resources, the VEP GRCh38 cache, Monarch KG, needLR backend files and MitoCarta, are not downloaded automatically by simply enabling a Snakemake rule. Check the configured paths before executing an analysis. Optional modules are controlled through the YAML settings; a configured path does not mean the analysis has completed.

## 3. Run the LRS workflow

Run these commands from the LRS workflow folder. Replace the Conda prefix if your environment location differs.

~~~bash
cd /path/to/SV-PIPELINE/snakemake_pipelines/lrs

# Preview planned jobs without executing them
snakemake -s Snakefile_LRS_update --configfile config_lrs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 --dry-run -p

# Run QC, callers, merged SVs and primary annotations
snakemake -s Snakefile_LRS_update --configfile config_lrs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 --rerun-incomplete -p

# Integrate evidence and generate downstream candidate tables
snakemake -s Snakefile_LRS_postprocess --configfile config_lrs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 -p

# Produce figures and plot-related data files
snakemake -s Snakefile_LRS_postprocess --configfile config_lrs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 all_thesis_plots
~~~

The postprocessing Snakefile normally selects samples with the required upstream files. Use the **postprocess_samples** configuration setting when you need to select specific ready samples. Snakemake can still schedule missing dependencies when a result is requested; always inspect a dry run before rebuilding earlier analyses.

## 4. What each LRS tool does

| Step / tool | Input | Role in the workflow |
| --- | --- | --- |
| **Dorado** | ONT BAM | Produces read-level sequencing-summary information. This workflow starts from existing BAMs; it does not perform the original basecalling. |
| **mosdepth** and BAM N50 script | Aligned BAM | Measure coverage, read lengths and sequencing quality. |
| **Clair3** | ONT BAM and reference | Calls genome-wide single-nucleotide variants (SNVs) and small insertions/deletions. |
| **WhatsHap** | Clair3 VCF and BAM | Phases small variants and generates a haplotagged BAM. |
| **Sniffles2, cuteSV, DELLY LR** | ONT BAM | Independently detect genome-wide structural variants, including deletions, insertions and rearrangements. |
| **Caller parser and evidence filter** | Caller VCFs | Put support fields into a common table and produce filtered caller VCFs before merging. |
| **Jasmine** | Filtered caller VCFs | Merges compatible calls into a complete patient-level SV VCF. A stricter multi-caller VCF is saved separately. |
| **needLR** | Dedicated Sniffles2-compatible query | Supplies long-read population-frequency context; it does not replace the Jasmine VCF. |
| **AnnotSV** | Merged SV VCF | Adds overlapping genes, SV annotation and available classification/population information. |
| **VEP** | Merged SV VCF | Adds transcript and sequence-consequence information. A panel-only view is made separately. |
| **Monarch KG / HPO** | Genes and phenotype resources | Links genes to phenotype evidence for research prioritization. |
| **Exomiser (offline)** | HON/patient HPO query | Provides a separate phenotype-based ranking for comparison. |
| **Straglr** | ONT BAM | Detects tandem-repeat variation. |
| **TLDR** | ONT BAM | Detects mobile-element insertion evidence. |
| **LongPhase** | ONT BAM and variants | Adds SNP/SV phasing information; its SNP and SV VCF outputs are distinct. |
| **modkit** | Modified-base ONT BAM | Summarizes CpG methylation. |
| **XCI analysis** | Chromosome-X phasing and methylation | Investigates X-inactivation skew where the required data and configuration exist. |
| **MitoCarta, gnomAD-SV and gene-depth scripts** | Annotated candidate tables | Add mitochondrial gene context, population comparison and local depth evidence. |
| **Samplot and plotting scripts** | BAM plus selected SVs and evidence tables | Visualize alignments, breakpoints, candidate loci and whole-genome summaries. |
| **Optional AlphaGenome module** | Small-variant VCF and configured AVI resource | Adds research-oriented predictions without changing variant calling or clinical classification. |

**Main distinction:** a call supported by several callers is not automatically pathogenic. Frequencies that are absent or cannot be matched remain unknown; they must not be interpreted as zero. Repeat, mobile-element, depth, methylation and phase results add context but are not independent experimental confirmation.

## 5. Main LRS outputs

These paths are inside the configured results folder, under each sample directory.

| Path under a sample | Meaning |
| --- | --- |
| **qc/**, **coverage/** | Dorado, BAM length/N50 and mosdepth QC. |
| **snp_clair3/** | Small-variant VCF. |
| **sv/sniffles2/**, **sv/cutesv/**, **sv/delly/** | Caller-specific SV VCFs. |
| **parsed/** | Standardized caller support and filtered evidence. |
| **sv/merged/SAMPLE_ID_merged_SV.vcf.gz** | Master Jasmine SV VCF (replace SAMPLE_ID with the sample). |
| **sv/merged/SAMPLE_ID_caller_support_summary.tsv** | Caller agreement and support. |
| **sv/needlr/** | Supplementary ONT population results. |
| **sv/annotsv/**, **sv/vep/** | Variant-to-gene and transcript annotation. |
| **gene_discovery/SAMPLE_ID_integrated_SV_gene_analysis.tsv** | Primary integrated SV–gene table. |
| **gene_discovery/SAMPLE_ID_allele_assessment.tsv** | Research allele and disease-mechanism assessment. |
| **gene_discovery/SAMPLE_ID_integrated_SV_gene_analysis.final.tsv.gz** | Later integrated analysis from postprocessing. |
| **gene_discovery/final/SAMPLE_ID_clinical_candidates.ranked.tsv** | Separate clinical candidate review order, not a diagnosis. |
| **phasing/**, **phasing_longphase/** | WhatsHap and LongPhase results. |
| **sv/straglr/**, **mei/tldr/**, **methylation/** | Tandem-repeat, mobile-element and methylation results. |
| **plots/** | Figures and plot-source tables when plotting is requested. |

Other final candidate, cohort, large-SV-depth and XCI tables are generated by the postprocessing and extension rules, depending on their inputs and enabled settings.

## 6. Shared Python and R scripts

The Snakefiles call most of these scripts automatically. Their names are listed here so you can find the code responsible for a particular result. For command-line options of a standalone Python script, run **python scripts/SCRIPT_NAME.py --help** from the repository root. Some modules are internal helpers and are intended to be imported, not run directly.

### Quality control and file checks

| Script | What it does |
| --- | --- |
| **scripts/summarize_dorado_qc.py** | Summarizes Dorado read metrics for ONT sequencing.. |
| **scripts/calc_bam_n50.py** | Calculates long-read length statistics and N50 directly from BAM alignments.. |
| **scripts/build_qc_summary.py** | Combines SRS mapping and coverage metrics into a QC table.. |
| **scripts/preflight_srs_wgs.py** | Checks SRS input files, genome reference and required resources before analysis.. |
| **scripts/audit_sv_outputs.py** | Checks that master SV records are accounted for in downstream results.. |
| **scripts/audit_vep_gene_coverage.py** | Compares gene coverage in VEP and AnnotSV annotations.. |

### SV calling, merging and evidence

| Script | What it does |
| --- | --- |
| **scripts/parse_sv_caller_vcf.py** | Converts caller-specific SV fields into a consistent evidence table.. |
| **scripts/filter_sv_evidence.py** | Applies configurable support, size and FILTER checks before merging.. |
| **scripts/summarize_sv_caller_support.py** | Reports which callers support each merged SV.. |
| **scripts/sanitize_jasmine_vcf.py** | Repairs missing VCF header definitions without changing variants.. |
| **scripts/combine_sv_cnv_evidence.py** | Joins SRS breakpoint SVs with CNVpytor read-depth evidence and retains passing depth-only CNVs.. |
| **scripts/add_supporting_evidence.py** | Adds GRIDSS, CNVpytor, MELT and repeat evidence to SRS candidate rows.. |
| **scripts/sv_evidence_common.py** | Shared parsers and helpers for conservative SV evidence handling.. |

### Annotation, genes and population resources

| Script | What it does |
| --- | --- |
| **scripts/build_integrated_sv_gene_tsv.py** | Combines master SVs, caller evidence, AnnotSV, population and gene annotations into one SV–gene table.. |
| **scripts/extract_annotsv_genes.py** | Extracts genes found in genome-wide AnnotSV results.. |
| **scripts/extract_nonpanel_genes.py** | Separates genes outside the optic-neuropathy panel.. |
| **scripts/filter_annotsv_panel.py** | Makes a panel-only view of complete AnnotSV results.. |
| **scripts/filter_vep_panel.py** | Makes a panel-only view of VEP results.. |
| **scripts/summarize_vep_transcripts.py** | Summarizes VEP transcript consequences, retaining the full annotation separately.. |
| **scripts/summarize_vep_sv_consequences.py** | Summarizes VEP SV consequences by variant and gene.. |
| **scripts/monarch_human_gene_phenotypes.py** | Gets gene–phenotype associations from a local Monarch Knowledge Graph.. |
| **scripts/annotate_gnomad_sv_exact.py** | Adds exact-site gnomAD-SV frequency matches where supported.. |
| **scripts/annotate_mitocarta_context.py** | Adds MitoCarta gene and pathway context.. |
| **scripts/annotate_straglr.py** | Links repeat-locus findings to overlapping genes.. |
| **scripts/annotsv_reconciliation.py** | Explains SV records omitted from AnnotSV output using the matching log.. |
| **scripts/audit_annotsv_evidence.py** | Reports annotation-source coverage and available AnnotSV database resources.. |
| **scripts/merge_integrated_annotation_deltas.py** | Combines supplemental annotation columns back into the main table.. |
| **scripts/summarize_secondary_mtdna_vcf.py** | Formats secondary ONT mtDNA small-variant findings.. |
| **scripts/summarize_mtdna_variants.py** | Formats SRS Mutserve mtDNA SNV and heteroplasmy results.. |

### Phenotype matching and candidate interpretation

| Script | What it does |
| --- | --- |
| **scripts/rank_sv_gene_candidates.py** | Prioritizes genes from genome-wide SV discovery using HON phenotype and gene-disease evidence.. |
| **scripts/build_hon_hpo_reference.py** | Builds an HON gene–disease–HPO reference from Monarch.. |
| **scripts/hpo_semantic_similarity.py** | Compares patient HPO terms with gene-associated phenotypes.. |
| **scripts/compare_exomiser_hon.py** | Runs or compares offline Exomiser phenotype matching against the HON method.. |
| **scripts/rank_sv_gene_events.py** | Ranks individual SV–gene relationships with the predicted variant effect.. |
| **scripts/build_candidate_tables.py** | Produces readable gene and SV–gene candidate tables.. |
| **scripts/assess_sv_alleles.py** | Organizes allele, disease-mechanism and missing-evidence information for research review.. |
| **scripts/build_diagnostic_candidates.py** | Creates SRS/LRS candidate review tables, including a separate clinical-review order.. |
| **scripts/finalize_candidate_ranking.py** | Combines later depth, phase, HPO and variant evidence for final review tables.. |
| **scripts/rank_mitochondrial_genes.py** | Prioritizes nuclear mitochondrial genes and separately reports mtDNA genes.. |
| **scripts/build_population_gene_effect_table.py** | Reports population frequencies alongside predicted SV–gene effects.. |
| **scripts/build_snv_sv_candidates.py** | Identifies small variants in genes that also contain an SV for manual second-hit review.. |
| **scripts/sv_gene_effects.py** | Describes how an SV overlaps, removes, duplicates or breaks near a gene.. |
| **scripts/ranking_common.py** | Shared logic and constants for ranking and evidence interpretation.. |

### Extra evidence and cohort comparisons

| Script | What it does |
| --- | --- |
| **scripts/summarize_large_sv_depth.py** | Summarizes read depth at and around large deletions and duplications.. |
| **scripts/build_large_sv_gene_context.py** | Lists genes intersecting large SVs for depth and locus review.. |
| **scripts/integrate_depth_into_candidates.py** | Adds depth observations to SV–gene candidate tables.. |
| **scripts/build_cohort_sv_recurrence.py** | Finds conservatively similar SVs across study samples.. |
| **scripts/filter_sv_candidates_for_review.py** | Creates separate size, frequency and supporting-evidence review views.. |
| **scripts/run_samplot_candidates.py** | Makes SRS read-alignment plots for selected SV candidates.. |
| **scripts/annotate_alphagenome_atlas.py** | Optional AlphaGenome offline variant-impact lookup and selected research plots.. |

### X-chromosome inactivation

| Script | What it does |
| --- | --- |
| **scripts/extract_xci_haplotags.py** | Extracts phased chromosome-X reads from a haplotagged BAM.. |
| **scripts/compare_xci_phasing.py** | Checks agreement between WhatsHap and LongPhase phase blocks on chromosome X.. |
| **scripts/xci_cluster_methylation.R** | Groups phased CpG methylation observations for XCI analysis.. |
| **scripts/calculate_xci_skew.py** | Estimates XCI skew and related confidence and sensitivity measures.. |

### Setup, rebuilds, checks and optional utilities

| Script | What it does |
| --- | --- |
| **scripts/check_annotsv_gene_resource.py** | Checks the prepared AnnotSV gene reference for invalid labels.. |
| **scripts/rebuild_annotsv_gene_resource.py** | Builds a corrected AnnotSV gene reference in a new folder.. |
| **scripts/install_annotsv_gene_resource.py** | Installs a validated corrected AnnotSV gene resource with backup.. |
| **scripts/check_run_cache.py** | Verifies saved inputs, settings and outputs before reusing expensive tool results.. |
| **scripts/rebuild_sv_results.py** | Recreates interpretation tables from saved results without rerunning SV callers.. |
| **scripts/run_cnvpytor_hg38.py** | Runs the GRCh38-specific CNVpytor compatibility wrapper.. |
| **scripts/run_melt_single.py** | Runs one MELT transposon family with a predictable output path.. |
| **scripts/benchmark_candidates.py** | Checks whether expected positive-control candidates were recovered.. |
| **scripts/benchmark_candidate_ranking.py** | Measures the rank of known variants or genes in a predefined truth set.. |
| **scripts/parse_vcf_to_tsv.py** | Standalone optional VCF-to-TSV inspection tool, not a Snakemake rule.. |

## 7. Plotting scripts

The **plots/run_thesis_plots.py** launcher uses the existing pipeline results to generate the figure collection. The individual figure scripts are grouped below. For outputs and options, see [plots/README.md](plots/README.md).

| Script | What it does |
| --- | --- |
| **plots/build_candidate_review.py** | Creates starter candidate-review resources. |
| **plots/build_gene_evidence_summary.py** | Summarizes gene-level evidence. |
| **plots/candidate_review.py** | Prepares candidate-review views. |
| **plots/integrate_candidate_context.py** | Adds supporting context to candidate views. |
| **plots/intersect_complementary_sv_evidence.py** | Matches repeat, MEI and other evidence by SV coordinates. |
| **plots/plot_caller_concordance.py** | Caller overlap plots. |
| **plots/plot_candidate_evidence_matrix.py** | Candidate evidence matrix. |
| **plots/plot_candidate_genes.py** | Prioritized gene plots. |
| **plots/plot_candidate_locus.py** | Per-candidate locus plots. |
| **plots/plot_cohort_sv_comparison.py** | Cross-sample SV comparison. |
| **plots/plot_downstream_sv_filtering.py** | Downstream review-filter plots. |
| **plots/plot_event_review_summary.py** | Event review summary. |
| **plots/plot_gene_hpo_heatmap.py** | Gene–HPO and phenotype-rank heatmaps. |
| **plots/plot_gene_multimodal_summary.py** | Combined gene-evidence figures. |
| **plots/plot_gene_sv_spectrum.py** | SV types associated with genes. |
| **plots/plot_genome_wide_sv_map.py** | Whole-genome SV positions and context. |
| **plots/plot_large_complex_candidates.py** | Large or complex SV overview. |
| **plots/plot_large_sv_overview.py** | Large-SV summaries. |
| **plots/plot_lrs_qc.py** | Long-read coverage figures. |
| **plots/plot_lrs_vs_srs.py** | LRS-versus-SRS benchmarking figures. |
| **plots/plot_mei.py** | Mobile-element insertion plots. |
| **plots/plot_methylation.py** | CpG methylation plots. |
| **plots/plot_mitocarta_sv_genes.py** | SVs affecting MitoCarta genes. |
| **plots/plot_mitochondrial_gene_ranking.py** | Mitochondrial gene rankings. |
| **plots/plot_needlr_population.py** | ONT population-frequency plots. |
| **plots/plot_panel_nonpanel_discovery.py** | Panel versus non-panel discovery plots. |
| **plots/plot_phasing_qc.py** | Phasing quality plots. |
| **plots/plot_straglr.py** | Tandem-repeat plots. |
| **plots/plot_sv_gene_associations.py** | SV–gene association plots. |
| **plots/plot_sv_gene_network.py** | SV-to-gene networks. |
| **plots/plot_sv_landscape.py** | SV distributions by class, size and location. |
| **plots/plot_utils.py** | Shared plotting helpers. |
| **plots/run_samplot_candidates.py** | ONT candidate read-alignment plots. |
| **plots/run_thesis_plots.py** | Launcher for the per-sample figure collection. |

The **r_plots/** folder contains complementary R figures:

| Script | What it does |
| --- | --- |
| **r_plots/plot_candidate_evidence_heatmap.R** | Candidate evidence heatmap. |
| **r_plots/plot_hpo_heatmap.R** | HPO association heatmap. |
| **r_plots/plot_large_sv_depth.R** | Large-SV depth plot. |
| **r_plots/plot_population_gene_effect.R** | Population and gene-effect plots. |
| **r_plots/plot_sequencing_qc.R** | Read-level and coverage QC figures. |
| **r_plots/plot_xci.R** | XCI methylation and skew figures. |
| **r_plots/theme_thesis.R** | Shared figure styling. |

## 8. How to read the results

Start with the master merged VCF and its caller-support summary, then check the integrated SV–gene table, gene/transcript effects, population evidence, HPO relevance, and the final candidate tables. Inspect the BAM or Samplot image for important breakpoints. The complete discovery variants are preserved even when an optional annotation is unavailable.

The gene-panel label and mitochondrial context help prioritize which findings to investigate; neither is evidence that an SV is disease-causing. An SV–gene overlap predicts a possible mechanism but does not demonstrate an effect on RNA or protein. Use the final TSV files to select variants for further analysis and experimental confirmation.

For short-read sequencing, follow the [SRS README](snakemake_pipelines/srs_wgs_pipeline/README.md). For the separate detailed checks and interpretation documents, consult **docs/**.
