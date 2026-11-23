# The confirmed truth set must have been curated BEFORE reviewing rankings.
# This optional rule runs after custom final ranking and Genomiser are complete.
if GENOMISER_ENABLED and config.get("genomiser_benchmark_truth_tsv"):
    rule benchmark_genomiser_ranking:
        input:
            truth=config["genomiser_benchmark_truth_tsv"],
            custom_genes=GENOMISER_CUSTOM_GENE,
            custom_sv=GENOMISER_CUSTOM_SV,
            genomiser_genes=PATH + "{sample}/genomiser/{sample}.genomiser.genes.tsv",
            genomiser_variants=PATH + "{sample}/genomiser/{sample}.genomiser.variants.tsv",
            script=SCRIPTS + "/benchmark_genomiser.py"
        output:
            details=PATH + "{sample}/benchmark/{sample}.genomiser_vs_custom.tsv",
            summary=PATH + "{sample}/benchmark/{sample}.genomiser_vs_custom.json"
        conda:
            CONDAENV + "genomiser.yaml"
        shell:
            """
            set -euo pipefail
            python {input.script:q} --sample {wildcards.sample:q}                 --truth {input.truth:q}                 --custom-genes {input.custom_genes:q} --custom-sv {input.custom_sv:q}                 --genomiser-genes {input.genomiser_genes:q}                 --genomiser-variants {input.genomiser_variants:q}                 --detail-output {output.details:q} --summary-output {output.summary:q}
            """
