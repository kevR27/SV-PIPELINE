# Non-destructive sensitivity of the FINAL SV-gene shortlist.
# This is not a truth benchmark or a replacement for clinical ranking.
if config.get("ranking_stability_enabled", True):
    rule ranking_stability:
        input:
            events=RANKING_STABILITY_SV,
            script=SCRIPTS + "/check_ranking_stability.py"
        output:
            summary=PATH + "{sample}/ranking_stability/stability_summary.tsv",
            genes=PATH + "{sample}/ranking_stability/stability_genes.tsv",
            manifest=PATH + "{sample}/ranking_stability/stability_manifest.json",
            pdf=PATH + "{sample}/ranking_stability/stability_summary.pdf",
            png=PATH + "{sample}/ranking_stability/stability_summary.png"
        params:
            outdir=PATH + "{sample}/ranking_stability",
            top_n=int(config.get("ranking_stability_top_n", 20))
        conda:
            CONDAENV + "plots.yaml"
        shell:
            """
            set -euo pipefail
            python {input.script:q} \\
                --input {input.events:q} \\
                --outdir {params.outdir:q} \\
                --top-n {params.top_n}
            test -s {output.manifest:q}
            """
