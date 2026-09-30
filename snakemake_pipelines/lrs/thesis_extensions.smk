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

VEP_CACHE_DIR_EXT = config.get("vep_cache_dir")
VEP_ASSEMBLY_EXT = config.get("vep_assembly", "GRCh38")

DORADO_BIN = config.get("dorado_bin", "dorado")
DORADO_SUMMARY_BAMS = config.get("dorado_summary_bams", {}) or {}
DORADO_SAMPLES = [sample for sample in POSTPROCESS_SAMPLES if sample in DORADO_SUMMARY_BAMS]
SNV_SV_SAMPLES = [
    sample
    for sample in POSTPROCESS_SAMPLES
    if os.path.exists(PATH + f"{sample}/phasing_longphase/{sample}.longphase.vcf.gz")
]


rule build_final_candidate_tables:
    input:
        genes=PATH + "{sample}/gene_discovery/{sample}_ranked_candidates.tsv",
        events=rules.rank_sv_gene_events.output.tsv,
        script=SCRIPTS + "/build_candidate_tables.py",
        effects=SCRIPTS + "/sv_gene_effects.py"
    output:
        genes=PATH + "{sample}/gene_discovery/final/{sample}_gene_candidates.tsv",
        sv=PATH + "{sample}/gene_discovery/final/{sample}_sv_gene_candidates.tsv"
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
        sv=rules.build_final_candidate_tables.output.sv,
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


rule r_candidate_evidence_heatmap:
    input:
        candidates=rules.build_final_candidate_tables.output.sv
    output:
        pdf=PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.pdf",
        png=PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.png",
        svg=PATH + "{sample}/plots_r/candidate_evidence/{sample}_candidate_evidence.svg"
    params:
        top_n=R_TOP_CANDIDATES
    conda:
        CONDAENV + "r_thesis_plots.yaml"
    script:
        "../../r_plots/plot_candidate_evidence_heatmap.R"


rule r_hpo_heatmap:
    input:
        phenotypes=PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv",
        genes=rules.build_final_candidate_tables.output.genes
    output:
        pdf=PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.pdf",
        png=PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.png",
        svg=PATH + "{sample}/plots_r/hpo/{sample}_hpo_heatmap.svg"
    params:
        top_genes=THESIS_TOP_GENES,
        top_hpo=R_TOP_HPO
    conda:
        CONDAENV + "r_thesis_plots.yaml"
    script:
        "../../r_plots/plot_hpo_heatmap.R"


rule r_large_sv_depth:
    input:
        bins=rules.summarize_large_sv_depth.output.bins,
        summary=rules.summarize_large_sv_depth.output.summary
    output:
        pdf=PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.pdf",
        png=PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.png",
        svg=PATH + "{sample}/plots_r/large_sv_depth/{sample}_large_sv_depth.svg"
    params:
        top_n=8
    conda:
        CONDAENV + "r_thesis_plots.yaml"
    script:
        "../../r_plots/plot_large_sv_depth.R"


def cohort_candidate_args():
    return " ".join(
        "--sample-input " + shlex.quote(
            sample + "=" + PATH + f"{sample}/gene_discovery/final/{sample}_sv_gene_candidates.tsv"
        )
        for sample in POSTPROCESS_SAMPLES
    )


rule cohort_sv_recurrence:
    input:
        tables=expand(
            PATH + "{sample}/gene_discovery/final/{sample}_sv_gene_candidates.tsv",
            sample=POSTPROCESS_SAMPLES,
        ),
        script=SCRIPTS + "/build_cohort_sv_recurrence.py"
    output:
        summary=PATH + "cohort_analysis/sv_recurrence.tsv",
        members=PATH + "cohort_analysis/sv_recurrence_members.tsv"
    params:
        sample_args=lambda wc: cohort_candidate_args()
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


rule dorado_read_summary:
    input:
        bam=lambda wc: DORADO_SUMMARY_BAMS[wc.sample]
    output:
        tsv=PATH + "{sample}/qc/dorado/{sample}.dorado_summary.tsv"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        {DORADO_BIN} summary {input.bam} > {output.tsv}
        test -s {output.tsv}
        """


rule summarize_dorado_qc:
    input:
        tsv=rules.dorado_read_summary.output.tsv,
        script=SCRIPTS + "/summarize_dorado_qc.py"
    output:
        tsv=PATH + "{sample}/qc/dorado/{sample}_dorado_qc.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script} \
            --input {input.tsv} \
            --sample {wildcards.sample} \
            --output {output.tsv}
        test -s {output.tsv}
        """


FINAL_THESIS_OUTPUTS = [
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
        PATH + "{sample}/gene_discovery/final/{sample}_snv_sv_candidates.tsv",
        sample=SNV_SV_SAMPLES,
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
        PATH + "{sample}/qc/dorado/{sample}_dorado_qc.tsv",
        sample=DORADO_SAMPLES,
    ),
]


rule final_thesis_analysis:
    input:
        FINAL_THESIS_OUTPUTS
