# Additional thesis/diagnostic-style analyses.
# This file is included by Snakefile_LRS_postprocess and does not replace
# any existing calling, annotation, ranking, or plotting rule.

R_PLOTS = os.path.abspath(os.path.join(workflow.basedir, "../../r_plots"))
FINAL_DIR = PATH + "{sample}/gene_discovery/final/"

DEPTH_BASE_BIN_BP = config.get("depth_base_bin_bp", 10000)
DEPTH_MIN_MAPQ = config.get("depth_min_mapq", 20)
LARGE_SV_DEPTH_MIN_SIZE = config.get("large_sv_depth_min_size", 100000)
R_TOP_CANDIDATES = config.get("r_plot_top_candidates", 30)
R_TOP_HPO = config.get("r_plot_top_hpo", 20)
HON_HPO_SEEDS = config.get(
    "hon_hpo_seed_terms",
    os.path.abspath(os.path.join(workflow.basedir, "../../reference/hon_hpo_seed_terms.tsv")),
)
HON_PANEL_FILE = config["candidate_genes_list"]
MONARCH_NODES_EXT = config["monarch_nodes"]
MONARCH_EDGES_EXT = config["monarch_edges"]
GENE_BED_EXT = config["gene_bed"]
PATIENT_HPO_FILES = config.get("patient_hpo_files", {}) or {}
GNOMAD_SV_VCF = config.get("gnomad_sv_vcf")


def patient_hpo_path(sample):
    value = PATIENT_HPO_FILES.get(str(sample))
    return str(value) if value else None

VEP_CACHE_DIR_EXT = config.get("vep_cache_dir")
VEP_ASSEMBLY_EXT = config.get("vep_assembly", "GRCh38")

# LongPhase is a mandatory LRS-update output. Use the same ready samples as the
# postprocess workflow rather than testing file existence at Snakefile parse
# time; parse-time existence checks can silently omit work in clean runs.
SNV_SV_SAMPLES = list(POSTPROCESS_SAMPLES)


rule build_hon_hpo_reference:
    input:
        panel=HON_PANEL_FILE,
        nodes=MONARCH_NODES_EXT,
        edges=MONARCH_EDGES_EXT,
        seeds=HON_HPO_SEEDS,
        pathways=MITOCARTA_PATHWAYS_FILE,
        script=SCRIPTS + "/build_hon_hpo_reference.py"
    output:
        tsv=PATH + "cohort_analysis/reference/hon_hpo_reference.tsv"
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        python {input.script} \
            --panel {input.panel} \
            --nodes {input.nodes} \
            --edges {input.edges} \
            --hon-seeds {input.seeds} \
            --mitopathways-gmx {input.pathways} \
            --output {output.tsv}
        test -s {output.tsv}
        """


def hpo_optional_input(wc):
    path = patient_hpo_path(wc.sample)
    return [path] if path else []


rule hpo_semantic_similarity:
    input:
        genes=PATH + "{sample}/gene_discovery/final/{sample}_gene_candidates.pre_depth.tsv",
        reference=rules.build_hon_hpo_reference.output.tsv,
        gene_phenotypes=PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv",
        edges=MONARCH_EDGES_EXT,
        patient=hpo_optional_input,
        script=SCRIPTS + "/hpo_semantic_similarity.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/final/{sample}_hpo_semantic_similarity.tsv"
    params:
        patient_arg=lambda wc: (
            f"--patient-hpo {shlex.quote(patient_hpo_path(wc.sample))}"
            if patient_hpo_path(wc.sample)
            else ""
        )
    conda:
        CONDAENV + "monarch.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --genes {input.genes} \
            --reference {input.reference} \
            --gene-phenotypes {input.gene_phenotypes} \
            --edges {input.edges} \
            {params.patient_arg} \
            --output {output.tsv}
        test -s {output.tsv}
        """


def gnomad_optional_input(wc):
    return [GNOMAD_SV_VCF] if GNOMAD_SV_VCF else []


# gnomAD-SV annotation is defined once in Snakefile_LRS_postprocess.\n# This extension consumes rules.annotate_gnomad_sv_exact.output.tsv.\n\nrule build_final_candidate_tables:
    input:
        genes=PATH + "{sample}/gene_discovery/{sample}_ranked_candidates.tsv",
        events=rules.annotate_gnomad_sv_exact.output.tsv,
        script=SCRIPTS + "/build_candidate_tables.py",
        effects=SCRIPTS + "/sv_gene_effects.py"
    output:
        genes=PATH + "{sample}/gene_discovery/final/{sample}_gene_candidates.pre_depth.tsv",
        sv=PATH + "{sample}/gene_discovery/final/{sample}_sv_gene_candidates.pre_depth.tsv"
    params:
        breakpoint_near=SV_GENE_BREAKPOINT_TOLERANCE
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.genes})
        python {input.script} \
            --gene-ranking {input.genes} \
            --sv-gene-events {input.events} \
            --gene-output {output.genes} \
            --sv-output {output.sv} \
            --breakpoint-near-bp {params.breakpoint_near}
        test -s {output.genes}
        test -s {output.sv}
        """


rule binned_read_depth:
    input:
        bam=lambda wc: BAMS[wc.sample]
    output:
        regions=PATH + "{sample}/coverage/binned/{sample}.10kb.regions.bed.gz",
        summary=PATH + "{sample}/coverage/binned/{sample}.10kb.mosdepth.summary.txt"
    params:
        prefix=PATH + "{sample}/coverage/binned/{sample}.10kb",
        bin_bp=DEPTH_BASE_BIN_BP,
        mapq=DEPTH_MIN_MAPQ
    threads:
        8
    conda:
        CONDAENV + "mosdepth.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.regions})
        mosdepth \
            --threads {threads} \
            --no-per-base \
            --by {params.bin_bp} \
            --mapq {params.mapq} \
            {params.prefix} \
            {input.bam}
        test -s {output.regions}
        test -s {output.summary}
        """


rule summarize_large_sv_depth:
    input:
        candidates=rules.build_final_candidate_tables.output.sv,
        coverage=rules.binned_read_depth.output.regions,
        script=SCRIPTS + "/summarize_large_sv_depth.py"
    output:
        summary=PATH + "{sample}/gene_discovery/final/{sample}_large_sv_depth.tsv",
        bins=PATH + "{sample}/gene_discovery/final/{sample}_large_sv_depth_bins.tsv"
    params:
        min_size=LARGE_SV_DEPTH_MIN_SIZE
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --candidates {input.candidates} \
            --coverage {input.coverage} \
            --summary-output {output.summary} \
            --bins-output {output.bins} \
            --min-size {params.min_size}
        test -e {output.summary}
        test -e {output.bins}
        """


rule build_large_sv_gene_context:
    input:
        candidates=rules.build_final_candidate_tables.output.sv,
        depth=rules.summarize_large_sv_depth.output.summary,
        gene_bed=GENE_BED_EXT,
        script=SCRIPTS + "/build_large_sv_gene_context.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/final/{sample}_large_sv_genes.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --candidates {input.candidates} \
            --depth-summary {input.depth} \
            --gene-bed {input.gene_bed} \
            --output {output.tsv}
        test -e {output.tsv}
        """


rule integrate_depth_into_final_candidates:
    input:
        genes=rules.build_final_candidate_tables.output.genes,
        sv=rules.build_final_candidate_tables.output.sv,
        depth=rules.summarize_large_sv_depth.output.summary,
        hpo=rules.hpo_semantic_similarity.output.tsv,
        script=SCRIPTS + "/integrate_depth_into_candidates.py"
    output:
        genes=PATH + "{sample}/gene_discovery/final/{sample}_gene_candidates.tsv",
        sv=PATH + "{sample}/gene_discovery/final/{sample}_sv_gene_candidates.tsv"
    params:
        min_size=LARGE_SV_DEPTH_MIN_SIZE
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --gene-candidates {input.genes} \
            --sv-candidates {input.sv} \
            --depth-summary {input.depth} \
            --hpo-similarity {input.hpo} \
            --gene-output {output.genes} \
            --sv-output {output.sv} \
            --min-size {params.min_size}
        test -s {output.genes}
        test -s {output.sv}
        """


rule vep_small_variants:
    input:
        vcf=PATH + "{sample}/phasing_longphase/{sample}.longphase.vcf.gz"
    output:
        txt=PATH + "{sample}/snp_clair3/{sample}.longphase.vep.txt"
    threads:
        4
    conda:
        CONDAENV + "vep.yaml"
    shell:
        """
        set -euo pipefail
        vep \
            --input_file {input.vcf} \
            --format vcf \
            --output_file {output.txt} \
            --cache \
            --offline \
            --dir_cache {VEP_CACHE_DIR_EXT} \
            --assembly {VEP_ASSEMBLY_EXT} \
            --canonical \
            --flag_pick \
            --symbol \
            --biotype \
            --variant_class \
            --force_overwrite \
            --fork {threads}
        test -s {output.txt}
        """


rule build_snv_sv_candidates:
    input:
        vep=rules.vep_small_variants.output.txt,
        phased=PATH + "{sample}/phasing_longphase/{sample}.longphase.vcf.gz",
        sv=rules.integrate_depth_into_final_candidates.output.sv,
        script=SCRIPTS + "/build_snv_sv_candidates.py"
    output:
        tsv=PATH + "{sample}/gene_discovery/final/{sample}_snv_sv_candidates.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --vep {input.vep} \
            --phased-vcf {input.phased} \
            --sv-candidates {input.sv} \
            --output {output.tsv}
        test -e {output.tsv}
        """


rule build_population_gene_effect:
    input:
        candidates=rules.integrate_depth_into_final_candidates.output.sv,
        script=SCRIPTS + "/build_population_gene_effect_table.py"
    output:
        table=PATH + "{sample}/gene_discovery/final/{sample}_population_gene_effect.tsv",
        genes=PATH + "{sample}/gene_discovery/final/{sample}_population_gene_summary.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --candidates {input.candidates} \
            --output {output.table} \
            --gene-summary {output.genes}
        test -s {output.table}
        test -s {output.genes}
        """


rule r_population_gene_effect:
    input:
        table=rules.build_population_gene_effect.output.table
    output:
        pdf=PATH + "{sample}/plots_r/population_gene_effect/{sample}_population_gene_effect.pdf",
        png=PATH + "{sample}/plots_r/population_gene_effect/{sample}_population_gene_effect.png",
        svg=PATH + "{sample}/plots_r/population_gene_effect/{sample}_population_gene_effect.svg"
    params:
        top_n=R_TOP_CANDIDATES
    conda:
        CONDAENV + "r_plot.yaml"
    script:
        "../../r_plots/plot_population_gene_effect.R"


rule r_candidate_evidence_heatmap:
    input:
        candidates=rules.integrate_depth_into_final_candidates.output.sv
    output:
        pdf=PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.pdf",
        png=PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.png",
        svg=PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.svg"
    params:
        top_n=R_TOP_CANDIDATES
    conda:
        CONDAENV + "r_plot.yaml"
    script:
        "../../r_plots/plot_candidate_evidence_heatmap.R"


rule r_hpo_heatmap:
    input:
        phenotypes=PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv",
        genes=rules.integrate_depth_into_final_candidates.output.genes
    output:
        pdf=PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.pdf",
        png=PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.png",
        svg=PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.svg"
    params:
        top_genes=THESIS_TOP_GENES,
        top_hpo=R_TOP_HPO
    conda:
        CONDAENV + "r_plot.yaml"
    script:
        "../../r_plots/plot_hpo_heatmap.R"


rule r_large_sv_depth:
    input:
        bins=rules.summarize_large_sv_depth.output.bins,
        summary=rules.summarize_large_sv_depth.output.summary,
        genes=rules.build_large_sv_gene_context.output.tsv
    output:
        pdf=PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.pdf",
        png=PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.png",
        svg=PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.svg"
    params:
        top_n=8
    conda:
        CONDAENV + "r_plot.yaml"
    script:
        "../../r_plots/plot_large_sv_depth.R"


def cohort_recurrence_args():
    return " ".join(
        "--sample-input " + shlex.quote(
            sample
            + "="
            + PATH
            + f"{sample}/gene_discovery/{sample}_integrated_SV_gene_with_mitocarta.tsv"
        )
        for sample in POSTPROCESS_SAMPLES
    )


rule cohort_sv_recurrence:
    input:
        tables=expand(
            PATH + "{sample}/gene_discovery/{sample}_integrated_SV_gene_with_mitocarta.tsv",
            sample=POSTPROCESS_SAMPLES,
        ),
        script=SCRIPTS + "/build_cohort_sv_recurrence.py"
    output:
        summary=PATH + "cohort_analysis/sv_recurrence.tsv",
        members=PATH + "cohort_analysis/sv_recurrence_members.tsv"
    params:
        sample_args=lambda wc: cohort_recurrence_args()
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {PATH}cohort_analysis
        python {input.script} \
            {params.sample_args} \
            --output {output.summary} \
            --members-output {output.members}
        test -s {output.summary}
        test -s {output.members}
        """


FINAL_THESIS_OUTPUTS = [
    rules.build_hon_hpo_reference.output.tsv,
    *expand(
        PATH + "{sample}/plots/.thesis_plots.done",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_gene_candidates.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_sv_gene_candidates.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_large_sv_depth.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_large_sv_genes.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_hpo_semantic_similarity.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_snv_sv_candidates.tsv",
        sample=SNV_SV_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_population_gene_effect.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/gene_discovery/final/{sample}_population_gene_summary.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/plots_r/population_gene_effect/{sample}_population_gene_effect.pdf",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.pdf",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.pdf",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.pdf",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/plots/13_breakpoint_evidence/" + SAMPLOT_RUN_NAME + "/.samplot.done",
        sample=POSTPROCESS_SAMPLES,
    ),
    PATH + "cohort_analysis/sv_recurrence.tsv",
    PATH + "cohort_analysis/sv_recurrence_members.tsv",
    *expand(
        PATH + "{sample}/qc/{sample}.dorado_qc.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
    *expand(
        PATH + "{sample}/qc/{sample}.n50.tsv",
        sample=POSTPROCESS_SAMPLES,
    ),
]


rule final_thesis_analysis:
    input:
        FINAL_THESIS_OUTPUTS
