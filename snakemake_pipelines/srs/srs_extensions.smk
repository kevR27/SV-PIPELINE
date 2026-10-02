# SRS-only evidence extensions.
# This file is intentionally not imported by the LRS workflows.

import os
import shlex


CNVPYTOR_RD_BIN = int(config.get("cnvpytor_rd_bin", 1000))
CNVPYTOR_BAF_BIN = int(config.get("cnvpytor_baf_bin", 10000))
CNVPYTOR_CHROMS = config.get(
    "cnvpytor_chromosomes",
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
)
CNVPYTOR_USE_1KG_MASK = bool(config.get("cnvpytor_use_1kg_mask", False))
CNVPYTOR_MAX_Q0 = float(config.get("cnvpytor_soft_max_q0", 0.5))
CNVPYTOR_MAX_PN = float(config.get("cnvpytor_soft_max_pn", 0.5))
CNVPYTOR_MIN_DG = float(config.get("cnvpytor_soft_min_gap_distance", 100000))

if CNVPYTOR_RD_BIN <= 0 or CNVPYTOR_RD_BIN % 100 != 0:
    raise ValueError("cnvpytor_rd_bin must be a positive multiple of 100")
if CNVPYTOR_BAF_BIN <= 0 or CNVPYTOR_BAF_BIN % 100 != 0:
    raise ValueError("cnvpytor_baf_bin must be a positive multiple of 100")

GRIDSS_BLACKLIST = config.get("gridss_blacklist_bed")
GRIDSS_JVM_HEAP_GB = int(config.get("gridss_jvm_heap_gb", 24))
GRIDSS_SKIP_SOFTCLIP = bool(config.get("gridss_skip_softclip_realignment", False))
GRIDSS_BWA_INDEXES = [REF + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]

MITY_REFERENCE = config.get("mity_reference", "hg38")
MITY_CONTIG = config.get("mity_contig", "chrM")
MITY_REPORT_MIN_VAF = float(config.get("mity_report_min_vaf", 0.01))


############################
# Additional Illumina QC   #
############################

rule samtools_wgs_qc:
    input:
        bam=lambda wc: BAMS[wc.sample],
        bai=lambda wc: bam_index_for(wc.sample)
    output:
        flagstat=PATH + "{sample}/qc/{sample}.flagstat.txt",
        stats=PATH + "{sample}/qc/{sample}.samtools.stats.txt",
        idxstats=PATH + "{sample}/qc/{sample}.idxstats.txt"
    threads: min(THREADS, 8)
    conda:
        CONDAENV + "srs_qc.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.flagstat})
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        samtools stats -@ {threads} {input.bam} > {output.stats}
        samtools idxstats {input.bam} > {output.idxstats}
        test -s {output.flagstat}
        test -s {output.stats}
        test -s {output.idxstats}
        """


rule summarize_srs_wgs_qc:
    input:
        flagstat=rules.samtools_wgs_qc.output.flagstat,
        stats=rules.samtools_wgs_qc.output.stats,
        mosdepth=rules.mosdepth_coverage.output.summary,
        script=SRS_SCRIPTS + "/summarize_srs_qc.py"
    output:
        tsv=PATH + "{sample}/qc/{sample}.srs_wgs_qc.tsv"
    conda:
        CONDAENV + "srs_qc.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} \
            --sample {wildcards.sample:q} \
            --flagstat {input.flagstat:q} \
            --stats {input.stats:q} \
            --mosdepth-summary {input.mosdepth:q} \
            --output {output.tsv:q}
        """


##############################################
# DeepVariant downstream WGS / panel views  #
##############################################

rule deepvariant_panel_subset:
    input:
        vcf=rules.deepvariant_snp.output.vcf,
        tbi=rules.deepvariant_snp.output.tbi,
        bed=PANEL_BED
    output:
        vcf=PATH + "{sample}/snp_deepvariant/{sample}.panel.vcf.gz",
        tbi=PATH + "{sample}/snp_deepvariant/{sample}.panel.vcf.gz.tbi"
    conda:
        CONDAENV + "srs_qc.yaml"
    shell:
        """
        set -euo pipefail
        bcftools view -R {input.bed:q} -Oz -o {output.vcf:q} {input.vcf:q}
        tabix -f -p vcf {output.vcf:q}
        """


##############################################
# CNVpytor: independent read-depth + BAF CNV #
##############################################

rule cnvpytor_snp_input:
    input:
        vcf=rules.deepvariant_snp.output.vcf,
        tbi=rules.deepvariant_snp.output.tbi
    output:
        vcf=PATH + "{sample}/cnv/cnvpytor/{sample}.biallelic_pass_snps.vcf.gz",
        tbi=PATH + "{sample}/cnv/cnvpytor/{sample}.biallelic_pass_snps.vcf.gz.tbi"
    conda:
        CONDAENV + "srs_cnvpytor.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.vcf})
        bcftools view -f PASS -m2 -M2 -v snps -Oz -o {output.vcf:q} {input.vcf:q}
        tabix -f -p vcf {output.vcf:q}
        """


rule cnvpytor_rd_baf:
    input:
        bam=lambda wc: BAMS[wc.sample],
        bai=lambda wc: bam_index_for(wc.sample),
        snps=rules.cnvpytor_snp_input.output.vcf,
        snps_tbi=rules.cnvpytor_snp_input.output.tbi
    output:
        root=PATH + "{sample}/cnv/cnvpytor/{sample}.pytor",
        raw_calls=PATH + "{sample}/cnv/cnvpytor/{sample}.calls.raw.tsv",
        genotypes=PATH + "{sample}/cnv/cnvpytor/{sample}.genotypes.raw.tsv"
    params:
        chroms=lambda wc: " ".join(shlex.quote(x) for x in CNVPYTOR_CHROMS),
        mask_flag=lambda wc: "" if CNVPYTOR_USE_1KG_MASK else "-nomask"
    conda:
        CONDAENV + "srs_cnvpytor.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.root})
        rm -f {output.root:q}

        cnvpytor -root {output.root:q} -rd {input.bam:q} -chrom {params.chroms}
        cnvpytor -root {output.root:q} -his {CNVPYTOR_RD_BIN}
        cnvpytor -root {output.root:q} -partition {CNVPYTOR_RD_BIN}
        cnvpytor -root {output.root:q} -call {CNVPYTOR_RD_BIN} > {output.raw_calls:q}

        cnvpytor -root {output.root:q} -snp {input.snps:q}
        if [ "{CNVPYTOR_USE_1KG_MASK}" = "True" ]; then
            cnvpytor -root {output.root:q} -mask_snps
        fi
        cnvpytor -root {output.root:q} -baf {CNVPYTOR_BAF_BIN} {params.mask_flag}

        awk 'NF >= 2 && ($1 == "deletion" || $1 == "duplication") {print $2}' {output.raw_calls:q} \
            | cnvpytor -root {output.root:q} -genotype {CNVPYTOR_BAF_BIN} -a {params.mask_flag} \
            > {output.genotypes:q}

        test -s {output.root}
        test -f {output.raw_calls}
        test -f {output.genotypes}
        """


rule parse_cnvpytor_calls:
    input:
        calls=rules.cnvpytor_rd_baf.output.raw_calls,
        genotypes=rules.cnvpytor_rd_baf.output.genotypes,
        script=SRS_SCRIPTS + "/parse_cnvpytor_calls.py"
    output:
        tsv=PATH + "{sample}/cnv/cnvpytor/{sample}.cnvpytor.calls.tsv"
    conda:
        CONDAENV + "srs_cnvpytor.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} \
            --calls {input.calls:q} \
            --genotypes {input.genotypes:q} \
            --sample {wildcards.sample:q} \
            --max-q0 {CNVPYTOR_MAX_Q0} \
            --max-pn {CNVPYTOR_MAX_PN} \
            --min-gap-distance {CNVPYTOR_MIN_DG} \
            --output {output.tsv:q}
        """


rule cnvpytor_to_vcf:
    input:
        tsv=rules.parse_cnvpytor_calls.output.tsv,
        ref=REF,
        fai=REF + ".fai",
        script=SRS_SCRIPTS + "/cnvpytor_calls_to_vcf.py"
    output:
        vcf=PATH + "{sample}/cnv/cnvpytor/{sample}.cnvpytor.vcf"
    conda:
        CONDAENV + "srs_cnvpytor.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} \
            --tsv {input.tsv:q} \
            --reference {input.ref:q} \
            --output {output.vcf:q}
        test -s {output.vcf}
        """


rule annotsv_cnvpytor:
    input:
        vcf=rules.cnvpytor_to_vcf.output.vcf,
        annotations=rules.annotsv_annotations_check.output.marker
    output:
        tsv=PATH + "{sample}/cnv/cnvpytor/{sample}.cnvpytor.annotsv.tsv"
    params:
        outdir=PATH + "{sample}/cnv/cnvpytor",
        hpo_flag=lambda wc: f"-hpo {shlex.quote(str(HPO_TERMS))}" if HPO_TERMS else ""
    conda:
        CONDAENV + "annotsv.yaml"
    shell:
        """
        set -euo pipefail
        AnnotSV \
            -SVinputFile {input.vcf:q} \
            -annotationsDir {ANNOTSV_ANNOTATIONS_DIR:q} \
            -outputDir {params.outdir:q} \
            -outputFile {wildcards.sample}.cnvpytor.annotsv.tsv \
            -genomeBuild GRCh38 \
            -candidateGenesFile {GENE_LIST:q} \
            -candidateGenesFiltering 0 \
            -tx ENSEMBL \
            -REreport 1 \
            -benignAF {ANNOTSV_BENIGN_AF} \
            {params.hpo_flag}
        test -f {output.tsv}
        """


#####################################################
# GRIDSS: independent breakpoint assembly evidence #
#####################################################

rule gridss_reference_check:
    input:
        ref=REF,
        fai=REF + ".fai"
    output:
        marker=PATH + "reference_checks/gridss_reference.ok"
    params:
        bwa_indexes=lambda wc: " ".join(shlex.quote(x) for x in GRIDSS_BWA_INDEXES)
    conda:
        CONDAENV + "srs_gridss.yaml"
    shell:
        """
        set -euo pipefail
        for f in {params.bwa_indexes}; do
            if [ ! -s "$f" ]; then
                echo "[GRIDSS] missing required BWA reference index: $f" >&2
                echo "[GRIDSS] create indexes with: bwa index {input.ref}" >&2
                exit 1
            fi
        done
        mkdir -p $(dirname {output.marker})
        touch {output.marker}
        """


rule gridss_breakpoints:
    input:
        bam=lambda wc: BAMS[wc.sample],
        bai=lambda wc: bam_index_for(wc.sample),
        ref=REF,
        fai=REF + ".fai",
        reference_ok=rules.gridss_reference_check.output.marker
    output:
        vcf=PATH + "{sample}/sv/gridss/{sample}.gridss.vcf.gz",
        tbi=PATH + "{sample}/sv/gridss/{sample}.gridss.vcf.gz.tbi",
        assembly=PATH + "{sample}/sv/gridss/{sample}.gridss.assembly.bam",
        assembly_bai=PATH + "{sample}/sv/gridss/{sample}.gridss.assembly.bam.bai"
    params:
        workdir=PATH + "{sample}/sv/gridss/work",
        blacklist=lambda wc: (
            "--blacklist " + shlex.quote(str(GRIDSS_BLACKLIST))
            if GRIDSS_BLACKLIST
            else ""
        ),
        skip_softclip=lambda wc: "--skipsoftcliprealignment" if GRIDSS_SKIP_SOFTCLIP else ""
    threads: THREADS
    conda:
        CONDAENV + "srs_gridss.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) {params.workdir}

        gridss \
            --reference {input.ref:q} \
            --output {output.vcf:q} \
            --assembly {output.assembly:q} \
            --threads {threads} \
            --workingdir {params.workdir:q} \
            --jvmheap {GRIDSS_JVM_HEAP_GB}g \
            {params.blacklist} \
            {params.skip_softclip} \
            {input.bam:q}

        bcftools index -f -t {output.vcf:q}
        samtools index -@ {threads} -f {output.assembly:q}
        test -s {output.vcf}
        test -s {output.tbi}
        test -s {output.assembly}
        test -s {output.assembly_bai}
        """


rule gridss_event_table:
    input:
        vcf=rules.gridss_breakpoints.output.vcf,
        tbi=rules.gridss_breakpoints.output.tbi,
        script=SRS_SCRIPTS + "/gridss_vcf_to_events.py"
    output:
        tsv=PATH + "{sample}/sv/gridss/{sample}.gridss.events.tsv"
    conda:
        CONDAENV + "srs_gridss.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} --vcf {input.vcf:q} --output {output.tsv:q}
        """


#############################################
# mity: dedicated mitochondrial WGS branch  #
#############################################

rule mity_call:
    input:
        bam=lambda wc: BAMS[wc.sample],
        bai=lambda wc: bam_index_for(wc.sample),
        fai=REF + ".fai"
    output:
        vcf=PATH + "{sample}/mtdna/mity/{sample}.mity.normalise.vcf.gz",
        tbi=PATH + "{sample}/mtdna/mity/{sample}.mity.normalise.vcf.gz.tbi"
    params:
        outdir=PATH + "{sample}/mtdna/mity"
    conda:
        CONDAENV + "srs_mity.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir:q}
        region=$(awk -v c="{MITY_CONTIG}" '$1 == c {print c ":1-" $2}' {input.fai:q})
        if [ -z "$region" ]; then
            echo "[mity] contig {MITY_CONTIG} not present in {input.fai}" >&2
            exit 1
        fi

        mity call \
            --reference {MITY_REFERENCE} \
            --prefix {wildcards.sample:q} \
            --output-dir {params.outdir:q} \
            --region "$region" \
            --normalise \
            {input.bam:q}

        test -s {output.vcf}
        test -s {output.tbi}
        """


rule mity_table:
    input:
        vcf=rules.mity_call.output.vcf,
        tbi=rules.mity_call.output.tbi,
        script=SRS_SCRIPTS + "/parse_mity_vcf.py"
    output:
        tsv=PATH + "{sample}/mtdna/mity/{sample}.mity.tsv"
    conda:
        CONDAENV + "srs_mity.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} \
            --vcf {input.vcf:q} \
            --sample {wildcards.sample:q} \
            --output {output.tsv:q}
        """


rule mity_report:
    input:
        vcf=rules.mity_call.output.vcf,
        tbi=rules.mity_call.output.tbi
    output:
        xlsx=PATH + "{sample}/mtdna/mity/{sample}.mity.report.xlsx"
    params:
        outdir=PATH + "{sample}/mtdna/mity"
    conda:
        CONDAENV + "srs_mity.yaml"
    shell:
        """
        set -euo pipefail
        mity report \
            --prefix {wildcards.sample:q} \
            --output-dir {params.outdir:q} \
            --min_vaf {MITY_REPORT_MIN_VAF} \
            {input.vcf:q}
        test -s {output.xlsx}
        """


#########################################################
# Final SRS evidence integration; master SVs preserved  #
#########################################################

rule augment_srs_evidence:
    input:
        integrated=rules.intersect_srs_orthogonal_evidence.output.tsv,
        cnvpytor=rules.parse_cnvpytor_calls.output.tsv,
        cnvpytor_annotsv=rules.annotsv_cnvpytor.output.tsv,
        gridss=rules.gridss_event_table.output.tsv,
        script=SRS_SCRIPTS + "/augment_srs_evidence.py"
    output:
        final=PATH + "{sample}/gene_discovery/{sample}_integrated_SRS_evidence.tsv",
        cnvpytor_only=PATH + "{sample}/gene_discovery/{sample}_cnvpytor_only_candidates.tsv",
        gridss_only=PATH + "{sample}/gene_discovery/{sample}_gridss_only_events.tsv"
    params:
        tolerance=ORTHOGONAL_BP_TOLERANCE
    conda:
        CONDAENV + "srs_cnvpytor.yaml"
    shell:
        """
        set -euo pipefail
        python {input.script:q} \
            --integrated {input.integrated:q} \
            --cnvpytor {input.cnvpytor:q} \
            --cnvpytor-annotsv {input.cnvpytor_annotsv:q} \
            --gridss-events {input.gridss:q} \
            --breakpoint-tolerance {params.tolerance} \
            --output {output.final:q} \
            --cnvpytor-only-output {output.cnvpytor_only:q} \
            --gridss-only-output {output.gridss_only:q}
        """
