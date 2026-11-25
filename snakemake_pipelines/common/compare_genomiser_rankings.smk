# Independent descriptive comparison for every analyzed sample (even unsolved).
rule compare_genomiser_rankings:
    input:
        custom_genes=GENOMISER_CUSTOM_GENE,
        genomiser_genes=PATH + "{sample}/genomiser/{sample}.genomiser.genes.tsv",
        script=SCRIPTS + "/compare_genomiser_rankings.py",
        helper=SCRIPTS + "/benchmark_genomiser.py"
    output:
        details=PATH + "{sample}/benchmark/{sample}.ranking_comparison.tsv",
        summary=PATH + "{sample}/benchmark/{sample}.ranking_comparison.json"
    conda:
        CONDAENV + "genomiser.yaml"
    shell:
        "python {input.script:q} --sample {wildcards.sample:q} --custom-genes {input.custom_genes:q} --genomiser-genes {input.genomiser_genes:q} --details {output.details:q} --summary {output.summary:q}"
