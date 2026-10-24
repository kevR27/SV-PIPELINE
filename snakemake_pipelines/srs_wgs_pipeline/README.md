# SRS whole-genome sequencing pipeline (Illumina)

This folder contains the **short-read sequencing (SRS)** workflow for analyzing Illumina whole-genome BAM files. It detects structural variants (SVs), small variants, copy-number changes and selected repeat/mobile-element variants, then connects these results to genes, phenotype information and supporting sequencing evidence.

SV and small-variant discovery are genome-wide. The hereditary optic-neuropathy (HON) panel and MitoCarta are used afterward to help organize the results. The mitochondrial genome is examined in a separate, secondary branch. This pipeline is for research prioritization, not automatic variant pathogenicity classification.

For long-read sequencing, see the [main LRS and shared-script README](../../README.md). For individual figures, see [plots/README.md](../../plots/README.md).

## 1. Files in this workflow

| File | How it is used |
| --- | --- |
| **Snakefile_SRS_WGS** | Main Snakemake workflow. Describes inputs, tools, outputs and how the analysis stages depend on one another. |
| **config_srs_wgs.yaml** | The sample BAM paths, output location, reference, annotation resources and run settings. This is the main file to edit. |
| **srs_interpretation.smk** | Downstream rules for SV–gene effects, depth, HPO matching, small-variant/SV pairs and final candidate tables. |
| **../../scripts/** | Shared Python/R helpers used by this Snakefile and the LRS workflow. |
| **../../envs/** | Conda definitions for the tools. |
| **../../plots/** and **../../r_plots/** | Figure and visualization scripts. |
| **../common/** | Shared rules for allele assessment, Exomiser, AnnotSV audits and optional AlphaGenome analysis. |

Each Snakemake **rule** represents one analysis step. **input** lists required files, **output** lists created files, **shell** runs the tool, and **conda** selects its environment. Snakemake follows file dependencies; the rule order on the page is not the job execution order.

## 2. What to prepare

The main input is a **coordinate-sorted, indexed paired-end Illumina WGS BAM**, aligned to the same GRCh38 reference specified in **config_srs_wgs.yaml**. The BAM, FASTA, annotation BEDs and other genomic resources must use matching reference coordinates.

Check these configuration fields before running:

| Setting | Meaning |
| --- | --- |
| **samples** | Sample name mapped to its BAM path. Quote numeric sample IDs. |
| **path** | Folder where the workflow writes the results. |
| **ref** | GRCh38 FASTA used for the BAM alignment. |
| **threads** | Maximum threads available for Snakemake rules. |
| **bed**, **candidate_genes_list**, **gene_bed** | HON panel and gene-location references, used during interpretation. |
| **vep_cache_dir**, **vep_cache_version**, **annotsv_annotations_dir** | Annotation data and matching versions. |
| **monarch_nodes**, **monarch_edges**, **hon_hpo_seed_terms** | Gene–phenotype/HPO resources. |
| **expansionhunter_catalog** | Repeat-locus catalog for ExpansionHunter. |
| **gridss_enabled**, **melt_enabled**, **mtdna_enabled**, **mitocarta_enabled** | Turn the indicated analysis branches on or off. |
| **deepvariant_use_gpu** | Switch between the configured CPU and GPU DeepVariant images. |
| **gnomad_sv_vcf**, **allele_patient_context**, **patient_hpo_files** | Optional population or patient-level interpretation input files. |

Conda environments come from **../../envs/**, but external databases, tool model files, VEP caches, MELT resources, GRIDSS reference preparation and other data may require manual setup. Setting an option to **true** does not install those resources. Check the matching environment definition and reference files first.

## 3. Run the pipeline

Start in this folder. Change the example Conda prefix if your environments are stored elsewhere.

~~~bash
cd /path/to/SV-PIPELINE/snakemake_pipelines/srs_wgs_pipeline

# Preview jobs and file dependencies
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 --dry-run -p

# Run the full analysis, including the downstream candidate tables
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 --rerun-incomplete -p all

# Run only the coverage/alignment QC target
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 quality_control

# Produce figures and the associated review/filter views
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 all_thesis_plots

# Generate downstream filtering tables without requesting all plots
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml --use-conda --conda-prefix /home/casadei7/snakemake_envs/envs/ --cores 32 further_filtering
~~~

The **all** target includes the main analysis and downstream interpretation. The **all_thesis_plots** target is separate. Requesting a specific output may still schedule missing upstream rules, so check a dry run before running only part of an older dataset.

## 4. Analysis steps and tools

The main workflow follows this sequence:

~~~text
Illumina WGS BAM
    |
    +-- Preflight, samtools and mosdepth -> alignment / coverage QC
    |
    +-- Manta + DELLY -> normalize and filter -> SURVIVOR
    |                                        |
    |                                        +-- merged breakpoint SV VCF
    |                                              |
    +-- GRIDSS ------------------------------------+-- additional breakpoint evidence
    |                                              |
    +-- CNVpytor -> copy-number evidence ----------+-- integrated SV VCF
    |                                                       |
    +-- DeepVariant -> VEP small variants                    +-- AnnotSV / VEP SV
    |                -> WhatsHap phasing                      +-- gene / phenotype context
    |                                                       +-- candidate review
    +-- ExpansionHunter / MELT / Mutserve (separate branches)
                                                            |
                                                  downstream interpretation
                                                            |
                                                  final tables and figures
~~~

### 4.1 Quality control

| Tool or script | What it does |
| --- | --- |
| **preflight_srs_wgs.py** | Checks that required input files, reference and resources can be found before expensive analyses run. |
| **samtools flagstat / stats** | Reports alignment and mapped-read statistics. |
| **mosdepth** | Calculates sequencing depth and coverage distributions, including depth bins used in later plots. |
| **build_qc_summary.py** | Combines the QC results into a readable sample summary. |

The QC step helps identify insufficient coverage, unusual mapping rates and reference/configuration problems. Its thresholds are review criteria, not evidence that all variant classes will be detectable.

### 4.2 Breakpoint structural variants

| Tool or script | What it does |
| --- | --- |
| **Manta** | Detects SVs from short-read pairing, split alignments and local assembly evidence. |
| **DELLY** | Detects SVs using paired-end and split-read information. |
| **parse_sv_caller_vcf.py** | Reads each caller's VCF and normalizes support and quality fields. |
| **filter_sv_evidence.py** | Applies configurable size, support and QC filters before merging. |
| **SURVIVOR** | Merges compatible Manta and DELLY calls into one breakpoint-SV VCF. |
| **summarize_sv_caller_support.py** | Reports whether a merged SV was found by Manta, DELLY or both. |
| **GRIDSS** | Provides separate breakpoint and assembly-based evidence, evaluated alongside the merged calls. |

The SURVIVOR master breakpoint VCF retains single-caller events. GRIDSS is not an additional SURVIVOR input. Agreement between callers from the same BAM is supporting computational evidence, not independent validation.

### 4.3 Small variants and phasing

| Tool or script | What it does |
| --- | --- |
| **DeepVariant** | Calls genome-wide SNVs and small insertions/deletions using the configured container image. |
| **VEP** | Annotates the full small-variant VCF, including noncoding consequences when supported by the configured cache. |
| **WhatsHap** | Phases compatible nearby small variants from the BAM and genotype calls. |
| **bcftools** | Selects a separate, focused nuclear mitochondrial review VCF without replacing the full VCF. |

The main VEP table retains genome-wide variant annotations. A small-variant result in a gene that also contains an SV does not automatically establish that the variants are on opposite chromosomes.

### 4.4 Copy-number variants and additional evidence

| Tool or script | What it does |
| --- | --- |
| **CNVpytor** | Detects deletions and duplications from read-depth changes using configured genome bins. |
| **run_cnvpytor_hg38.py** | Runs the repository's GRCh38-specific CNVpytor wrapper. |
| **combine_sv_cnv_evidence.py** | Adds CNVpytor depth information to breakpoint SVs and retains passing CNVpytor-only CNVs. |
| **add_supporting_evidence.py** | Adds depth, GRIDSS, repeat-locus and optional MELT evidence to candidate tables. |

The **integrated SV VCF** is built from the SURVIVOR breakpoint calls and passing copy-number evidence. CNVpytor-only events may need extra work to define exact breakpoints and genotype.

### 4.5 Annotation and phenotype interpretation

| Tool or script | What it does |
| --- | --- |
| **AnnotSV** | Identifies overlapping genes and supplies relevant structural-variant annotation and reported classification fields. |
| **VEP SV annotation** | Describes transcript-level consequences and affected features. |
| **extract_annotsv_genes.py**, **extract_nonpanel_genes.py** | Lists genes affected by SVs and distinguishes panel/non-panel genes. |
| **filter_annotsv_panel.py**, **filter_vep_panel.py** | Creates panel-specific views without limiting genome-wide discovery. |
| **monarch_human_gene_phenotypes.py** | Retrieves gene–HPO associations from locally downloaded Monarch data. |
| **rank_sv_gene_candidates.py** | Prioritizes genes for the HON research question. |
| **compare_exomiser_hon.py** | Uses offline Exomiser phenotype ranking as a separate comparison. |
| **build_integrated_sv_gene_tsv.py** | Creates the central SV–gene evidence table from the complete integrated SV VCF and its annotations. |
| **assess_sv_alleles.py** | Reports allele and disease-mechanism hypotheses, supporting evidence and unknown information. |
| **build_diagnostic_candidates.py** | Prepares candidate review tables and a separate clinical review order. |
| **annotate_mitocarta_context.py** | Adds mitochondrial gene membership and pathway information from MitoCarta. |
| **rank_mitochondrial_genes.py** | Organizes nuclear-encoded mitochondrial candidates and secondary mtDNA context. |
| **annotate_gnomad_sv_exact.py** | Looks for conservative exact-site population matches when gnomAD-SV is configured. |

The HON panel and MitoCarta are annotation layers, not evidence that an observed variant is disease-causing. Missing population information is kept separate from a measured low allele frequency.

### 4.6 Other variant types

| Tool or script | What it does |
| --- | --- |
| **ExpansionHunter** | Genotypes the repeat loci described in the configured variant catalog. |
| **MELT** | Identifies selected mobile-element insertion families when enabled and prepared. |
| **run_melt_single.py** | Runs an individual MELT family and gives its output a predictable filename. |
| **Mutserve** | Calls mtDNA SNVs and estimates heteroplasmy for the secondary mitochondrial branch. |
| **summarize_mtdna_variants.py** | Converts Mutserve results into a more readable TSV. |

ExpansionHunter is catalog-based and is not a complete genome-wide repeat-discovery workflow. MELT and Mutserve do not replace the primary nuclear SV analysis.

## 5. Downstream analysis: the interpretation module

**srs_interpretation.smk** is included by the main Snakefile. It does not repeat variant calling. Instead, it uses completed VCFs, annotations, gene tables and coverage results to create a more detailed review.

| Script | Output / purpose |
| --- | --- |
| **rank_sv_gene_events.py** | Lists individual SV–gene relationships and possible variant effects. |
| **build_candidate_tables.py** | Makes separate readable gene and SV–gene candidate tables. |
| **summarize_large_sv_depth.py** | Checks local depth around large deletions and duplications. |
| **build_large_sv_gene_context.py** | Lists all genes in large SV intervals. |
| **build_hon_hpo_reference.py**, **hpo_semantic_similarity.py** | Build phenotype references and compare patient HPO terms with candidate genes when supplied. |
| **integrate_depth_into_candidates.py** | Adds large-SV depth observations to candidate results. |
| **build_snv_sv_candidates.py** | Flags small variants and SVs affecting the same gene for second-hit investigation. |
| **finalize_candidate_ranking.py** | Builds final review tables from available interpretation evidence. |

These rules use the shared Python functions in **../../scripts/ranking_common.py**, **../../scripts/sv_gene_effects.py** and **../../scripts/sv_evidence_common.py**. The final ranking does not replace the original calls, AnnotSV classification or separate Exomiser results.

For the complete shared-script inventory, including setup/rebuild utilities and the LRS-only XCI scripts, see [the main README](../../README.md#6-shared-python-and-r-scripts).

## 6. Where to find the results

Every sample gets its own directory inside the **path** configured in the YAML file. Replace SAMPLE_ID below with the actual sample name.

| Folder or file under SAMPLE_ID/ | What it contains |
| --- | --- |
| **qc/SAMPLE_ID.wgs_qc_summary.tsv** | QC summary. |
| **coverage/** | mosdepth coverage reports and depth bins. |
| **sv/manta/**, **sv/delly/**, **sv/gridss/** | Caller VCFs and GRIDSS evidence, if enabled. |
| **sv/merged/SAMPLE_ID_merged_SV.vcf.gz** | Manta/DELLY SURVIVOR breakpoint-SV VCF. |
| **cnv/cnvpytor/** | Copy-number calls, statistics and CNVpytor files. |
| **sv/integrated/SAMPLE_ID_integrated_SV.vcf.gz** | Main combined breakpoint-SV and passing CNV discovery VCF. |
| **snp_deepvariant/SAMPLE_ID.vcf.gz** | Genome-wide DeepVariant small variants. |
| **snp_deepvariant/SAMPLE_ID.genomewide.vep.tsv.gz** | Genome-wide VEP small-variant annotations. |
| **phasing/** | WhatsHap small-variant phasing. |
| **sv/annotsv/** and **sv/vep/** | SV gene and transcript annotations. |
| **gene_discovery/SAMPLE_ID_integrated_SV_gene_analysis.tsv** | Integrated SV–gene evidence. |
| **gene_discovery/SAMPLE_ID_allele_assessment.tsv** | Research allele-level assessment. |
| **diagnostic_review/SAMPLE_ID.diagnostic_candidates.tsv** | Main initial candidate review table. |
| **diagnostic_review/SAMPLE_ID.clinical_candidates.tsv** | Separate clinical candidate review order. |
| **gene_discovery/final/** | Later SV–gene events, large-SV depth, SNV–SV comparisons and final rankings. |
| **str/expansionhunter/**, **mei/melt/**, **mtDNA/** | Separate repeat, mobile-element and mitochondrial analysis results, when enabled. |
| **diagnostic_review/samplot/** | Read-alignment images for selected candidates. |
| **downstream_filtering/** and **plots/** | Additional review views and generated figures. |

You can request **all_thesis_plots** to make figures after the analysis. The plotting launcher is **../../plots/run_thesis_plots.py**, with individual plotting scripts described in the [main README](../../README.md#7-plotting-scripts) and [plotting guide](../../plots/README.md). SRS review/filter plots use **../../plots/plot_downstream_sv_filtering.py** and **../../scripts/filter_sv_candidates_for_review.py**.

## 7. Interpreting the findings

Begin with the QC summary and the integrated SV VCF, then inspect the SV–gene table, variant effect, caller support, read depth, population evidence and phenotype associations. For high-priority events, inspect the alignment with Samplot or a BAM viewer, especially when a deletion/duplication may represent a complex rearrangement.

The review outputs distinguish **research ranking**, **reported AnnotSV classification** and **clinical review order**. They are not equivalent. A variant may overlap a relevant gene while still requiring breakpoint confirmation, segregation analysis, transcript studies or other experimental follow-up. The workflow does not infer inheritance, variant phase or pathogenicity from missing evidence.
