# Optional X-chromosome inactivation (XCI) analysis.
#
# This layer is modeled on the read-level logic used by SkewX/Gocuk et al.:
#   native-DNA CpG methylation -> two epiallele clusters (Xa/Xi)
#   + WhatsHap HP/PS read tags -> block-wise skew
#   + folded-binomial maximum-likelihood global skew.
#
# LongPhase is used as an independent phase-consistency QC layer. The analysis
# is opt-in because XCI skew is only interpretable for appropriate XX samples
# with native-DNA modified-base calls and sufficient chrX coverage.

XCI_SAMPLES = [str(x) for x in (config.get("xci_samples") or [])]
XCI_UNKNOWN = sorted(set(XCI_SAMPLES) - set(str(x) for x in SAMPLES))
if XCI_UNKNOWN:
    raise ValueError(
        "xci_samples contains sample IDs not present in config samples: "
        + ", ".join(XCI_UNKNOWN)
    )

XCI_CHROM = str(config.get("xci_chrom", "chrX"))
XCI_CPG_ISLANDS_BED = config.get("xci_cpg_islands_bed")
XCI_XIST_PROMOTER_BED = config.get("xci_xist_promoter_bed")
XCI_MIN_CLUSTER_READS = int(config.get("xci_min_cluster_reads", 5))
XCI_MIN_CHRX_COVERAGE = float(config.get("xci_min_chrx_coverage", 15))
XCI_MIN_BLOCK_READS = int(config.get("xci_min_block_reads", 5))
XCI_NEUTRAL_THRESHOLD = float(config.get("xci_neutral_threshold", 0.40))
XCI_HIGH_SKEW_THRESHOLD = float(config.get("xci_high_skew_threshold", 0.20))
XCI_METHYLATION_BIN_BP = int(config.get("xci_methylation_bin_bp", 5000000))
XCI_THREADS = int(config.get("xci_threads", min(int(config.get("threads", 8)), 8)))
XCI_MODKIT = config.get("modkitenv", "modkit")
XCI_R_PLOTS = os.path.join(REPO_ROOT, "r_plots")
XCI_MODKIT_IO_THREADS = int(config.get("modkit_io_threads", XCI_THREADS))
XCI_MODKIT_SAMPLING_THREADS = int(
    config.get("modkit_sampling_threads", XCI_THREADS)
)
XCI_MODKIT_BGZF_THREADS = int(config.get("modkit_bgzf_threads", XCI_THREADS))

if XCI_SAMPLES and not XCI_CPG_ISLANDS_BED:
    raise ValueError(
        "xci_cpg_islands_bed must be configured when xci_samples is non-empty. "
        "Use a GRCh38 chrX CpG-island BED because the active LRS reference is hg38."
    )

XCI_FINAL_OUTPUTS = [
    PATH + f"{sample}/xci/{sample}_xci_summary.tsv"
    for sample in XCI_SAMPLES
]
XCI_PLOT_DONE = [
    PATH + f"{sample}/plots/15_x_inactivation/.xci_plots.done"
    for sample in XCI_SAMPLES
]


def xci_xist_input(wc):
    return [XCI_XIST_PROMOTER_BED] if XCI_XIST_PROMOTER_BED else []


def xci_xist_arg(wc):
    return (
        shlex.quote(str(XCI_XIST_PROMOTER_BED))
        if XCI_XIST_PROMOTER_BED
        else "NONE"
    )


rule xci_extract_haplotags:
    input:
        bam=PATH + "{sample}/phasing/{sample}.phased.bam",
        bai=PATH + "{sample}/phasing/{sample}.phased.bam.bai"
    output:
        reads=PATH + "{sample}/xci/{sample}_chrX_whatshap_reads.tsv.gz",
        summary=PATH + "{sample}/xci/{sample}_chrX_whatshap_haplotag_summary.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {PATH}{wildcards.sample}/xci

        python {SCRIPTS}/extract_xci_haplotags.py \
            --bam {input.bam} \
            --chrom {XCI_CHROM} \
            --output {output.reads} \
            --summary-output {output.summary}

        test -s {output.reads}
        test -s {output.summary}
        """


rule xci_compare_phasing:
    input:
        whatshap=PATH + "{sample}/phasing/{sample}.phased.vcf.gz",
        whatshap_index=PATH + "{sample}/phasing/{sample}.phased.vcf.gz.tbi",
        longphase=PATH + "{sample}/phasing_longphase/{sample}.longphase.vcf.gz",
        longphase_index=PATH + "{sample}/phasing_longphase/{sample}.longphase.vcf.gz.tbi"
    output:
        blocks=PATH + "{sample}/xci/{sample}_whatshap_longphase_phase_blocks.tsv",
        summary=PATH + "{sample}/xci/{sample}_whatshap_longphase_phase_summary.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {PATH}{wildcards.sample}/xci

        python {SCRIPTS}/compare_xci_phasing.py \
            --whatshap-vcf {input.whatshap} \
            --longphase-vcf {input.longphase} \
            --chrom {XCI_CHROM} \
            --blocks-output {output.blocks} \
            --summary-output {output.summary}

        test -s {output.blocks}
        test -s {output.summary}
        """


rule xci_modkit_phased_chrX:
    input:
        bam=PATH + "{sample}/phasing/{sample}.phased.bam",
        bai=PATH + "{sample}/phasing/{sample}.phased.bam.bai",
        ref=REF,
        haplotag_summary=rules.xci_extract_haplotags.output.summary
    output:
        hp1=PATH + "{sample}/xci/methylation/{sample}.chrX.hp1.5mC.bedmethyl.gz",
        hp2=PATH + "{sample}/xci/methylation/{sample}.chrX.hp2.5mC.bedmethyl.gz",
        combined=PATH + "{sample}/xci/methylation/{sample}.chrX.combined.5mC.bedmethyl.gz"
    log:
        PATH + "{sample}/logs/xci/{sample}.modkit_phased.log"
    threads:
        XCI_THREADS
    conda:
        CONDAENV + "modkit.yaml"
    shell:
        """
        set -euo pipefail

        outdir={PATH}{wildcards.sample}/xci/methylation
        tmpdir="$outdir/.{wildcards.sample}.modkit_phased_tmp"

        mkdir -p "$outdir"
        mkdir -p $(dirname {log})
        rm -rf "$tmpdir"
        mkdir -p "$tmpdir"

        # Check that WhatsHap produced usable chrX haplotags before asking
        # modkit for phased methylation. XCI needs reads from both haplotypes
        # and phase-set information; otherwise hp1/hp2 output will be empty.
        hp1_reads=$(awk -F'\\t' '
            NR==1 {for(i=1;i<=NF;i++) if($i=="HP1_reads") c=i; next}
            NR==2 && c {print $c}
        ' {input.haplotag_summary})
        hp2_reads=$(awk -F'\\t' '
            NR==1 {for(i=1;i<=NF;i++) if($i=="HP2_reads") c=i; next}
            NR==2 && c {print $c}
        ' {input.haplotag_summary})
        ps_reads=$(awk -F'\\t' '
            NR==1 {for(i=1;i<=NF;i++) if($i=="reads_with_PS") c=i; next}
            NR==2 && c {print $c}
        ' {input.haplotag_summary})

        if [[ -z "$hp1_reads" || -z "$hp2_reads" || -z "$ps_reads" ]]; then
            echo "[ERROR] Could not read HP1_reads, HP2_reads or reads_with_PS from {input.haplotag_summary}" >&2
            cat {input.haplotag_summary} >&2
            exit 1
        fi

        if (( hp1_reads == 0 || hp2_reads == 0 || ps_reads == 0 )); then
            echo "[ERROR] XCI cannot use {wildcards.sample}: chrX WhatsHap haplotags are insufficient." >&2
            echo "        HP1_reads=$hp1_reads HP2_reads=$hp2_reads reads_with_PS=$ps_reads" >&2
            echo "        Check that the sample is appropriate for XCI and that chrX has enough heterozygous variants for phasing." >&2
            exit 1
        fi

        echo "[XCI] chrX WhatsHap QC: HP1=$hp1_reads HP2=$hp2_reads reads_with_PS=$ps_reads"

        # Primary XCI methylation analysis uses 5mCG only, matching the
        # methylation mark used in the published nanopore XCI workflow.
        # --phased partitions the HP-tagged WhatsHap modBAM into HP1/HP2.
        {XCI_MODKIT} pileup \
            {input.bam} \
            "$tmpdir" \
            --reference {input.ref} \
            --region {XCI_CHROM} \
            --cpg \
            --modified-bases 5mC \
            --combine-strands \
            --phased \
            --bgzf \
            --threads {threads} \
            --io-threads {XCI_MODKIT_IO_THREADS} \
            --sampling-threads {XCI_MODKIT_SAMPLING_THREADS} \
            --bgzf-threads {XCI_MODKIT_BGZF_THREADS} \
            --log-filepath {log}

        if [[ ! -s "$tmpdir/hp1.bedmethyl.gz" || ! -s "$tmpdir/hp2.bedmethyl.gz" || ! -s "$tmpdir/combined.bedmethyl.gz" ]]; then
            echo "[ERROR] modkit finished but one or more phased bedMethyl outputs are missing or empty." >&2
            echo "        Expected: hp1.bedmethyl.gz, hp2.bedmethyl.gz, combined.bedmethyl.gz" >&2
            ls -lah "$tmpdir" >&2 || true
            echo "        WhatsHap counts before modkit: HP1=$hp1_reads HP2=$hp2_reads PS=$ps_reads" >&2
            exit 1
        fi

        mv "$tmpdir/hp1.bedmethyl.gz" {output.hp1}
        mv "$tmpdir/hp2.bedmethyl.gz" {output.hp2}
        mv "$tmpdir/combined.bedmethyl.gz" {output.combined}
        rm -rf "$tmpdir"

        test -s {output.hp1}
        test -s {output.hp2}
        test -s {output.combined}
        """


rule xci_cluster_methylation:
    input:
        bam=PATH + "{sample}/phasing/{sample}.phased.bam",
        cpg=lambda wc: XCI_CPG_ISLANDS_BED,
        haplotags=rules.xci_extract_haplotags.output.reads,
        xist=xci_xist_input
    output:
        clustered=PATH + "{sample}/xci/{sample}_chrX_clustered_methylation_reads.tsv.gz",
        raw_blocks=PATH + "{sample}/xci/{sample}_chrX_block_skew.raw.tsv.gz"
    params:
        xist=xci_xist_arg
    threads:
        XCI_THREADS
    conda:
        CONDAENV + "xci.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {PATH}{wildcards.sample}/xci

        Rscript {SCRIPTS}/xci_cluster_methylation.R \
            {wildcards.sample} \
            {input.bam} \
            {input.cpg} \
            {input.haplotags} \
            {output.clustered} \
            {output.raw_blocks} \
            {XCI_MIN_CLUSTER_READS} \
            {threads} \
            {params.xist}

        test -s {output.clustered}
        test -s {output.raw_blocks}
        """


rule xci_calculate_skew:
    input:
        blocks=rules.xci_cluster_methylation.output.raw_blocks,
        coverage=PATH + "{sample}/coverage/{sample}.mosdepth.summary.txt",
        haplotags=rules.xci_extract_haplotags.output.summary,
        phase=rules.xci_compare_phasing.output.summary
    output:
        blocks=PATH + "{sample}/xci/{sample}_xci_blocks.tsv",
        summary=PATH + "{sample}/xci/{sample}_xci_summary.tsv"
    conda:
        CONDAENV + "plots.yaml"
    shell:
        """
        set -euo pipefail

        python {SCRIPTS}/calculate_xci_skew.py \
            --block-skew {input.blocks} \
            --mosdepth-summary {input.coverage} \
            --haplotag-summary {input.haplotags} \
            --phase-summary {input.phase} \
            --chrom {XCI_CHROM} \
            --min-chrx-coverage {XCI_MIN_CHRX_COVERAGE} \
            --min-block-reads {XCI_MIN_BLOCK_READS} \
            --neutral-threshold {XCI_NEUTRAL_THRESHOLD} \
            --high-skew-threshold {XCI_HIGH_SKEW_THRESHOLD} \
            --blocks-output {output.blocks} \
            --summary-output {output.summary}

        test -s {output.blocks}
        test -s {output.summary}
        """


rule xci_plot_analysis:
    input:
        blocks=rules.xci_calculate_skew.output.blocks,
        summary=rules.xci_calculate_skew.output.summary,
        phase_blocks=rules.xci_compare_phasing.output.blocks,
        phase_summary=rules.xci_compare_phasing.output.summary,
        hp1=rules.xci_modkit_phased_chrX.output.hp1,
        hp2=rules.xci_modkit_phased_chrX.output.hp2,
        script=XCI_R_PLOTS + "/plot_xci.R"
    output:
        done=PATH + "{sample}/plots/15_x_inactivation/.xci_plots.done"
    params:
        outdir=lambda wc: (
            PATH + f"{wc.sample}/plots/15_x_inactivation"
        )
    conda:
        CONDAENV + "r_plot.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        Rscript {input.script} \
            {input.blocks} \
            {input.summary} \
            {input.phase_blocks} \
            {input.phase_summary} \
            {input.hp1} \
            {input.hp2} \
            {params.outdir} \
            {wildcards.sample} \
            {XCI_METHYLATION_BIN_BP}

        test -s {params.outdir}/{wildcards.sample}_xci_block_skew_distribution.pdf
        test -s {params.outdir}/{wildcards.sample}_xci_chrX_block_skew.pdf
        test -s {params.outdir}/{wildcards.sample}_xci_orientation_log_odds.pdf
        test -s {params.outdir}/{wildcards.sample}_xci_haplotype_methylation.pdf
        test -s {params.outdir}/{wildcards.sample}_xci_phase_concordance.pdf
        test -s {params.outdir}/{wildcards.sample}_xci_plot_manifest.tsv

        touch {output.done}
        """


rule xci_analysis:
    input:
        XCI_FINAL_OUTPUTS,
        XCI_PLOT_DONE
