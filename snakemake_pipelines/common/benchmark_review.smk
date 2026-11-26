# Per-patient benchmarking index. Distinguish confirmed truth, rank agreement
# and descriptive filtering sensitivity; never combine them into one accuracy metric.
rule benchmark_review:
    input:
        comparison=PATH + "{sample}/benchmark/{sample}.ranking_comparison.json",
        stability=lambda wc: (
            [PATH + f"{wc.sample}/ranking_stability/stability_summary.tsv"]
            if config.get("ranking_stability_enabled", True) else []
        ),
        truth=lambda wc: (
            [PATH + f"{wc.sample}/benchmark/{wc.sample}.genomiser_vs_custom.json"]
            if wc.sample in GENOMISER_BENCHMARK_SAMPLES else []
        ),
        script=SCRIPTS + "/summarize_benchmarking_review.py"
    output:
        table=PATH + "{sample}/benchmark/{sample}.benchmark_review.tsv",
        manifest=PATH + "{sample}/benchmark/{sample}.benchmark_review.json"
    params:
        stability_arg=lambda wc: (
            "--stability " + shlex.quote(PATH + f"{wc.sample}/ranking_stability/stability_summary.tsv")
            if config.get("ranking_stability_enabled", True) else ""
        ),
        truth_arg=lambda wc: (
            "--truth " + shlex.quote(PATH + f"{wc.sample}/benchmark/{wc.sample}.genomiser_vs_custom.json")
            if wc.sample in GENOMISER_BENCHMARK_SAMPLES else ""
        )
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --sample {wildcards.sample:q} --comparison {input.comparison:q} {params.stability_arg} {params.truth_arg} --output {output.table:q} --manifest {output.manifest:q}
        test -s {output.table:q}
        test -s {output.manifest:q}
        """

rule benchmarking_review_all:
    input:
        expand(PATH + "{sample}/benchmark/{sample}.benchmark_review.tsv",
               sample=BENCHMARK_REVIEW_SAMPLES)
