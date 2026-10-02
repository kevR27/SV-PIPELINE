# SRS WGS v2 — Illumina whole-genome structural-variant workflow

This folder is an isolated successor to the original SRS workflow. It does not
modify the LRS workflow or the original SRS Snakefile.

## Design principle

The workflow separates **variant/event discovery** from **independent evidence
modalities**. Manta and DELLY define the primary breakpoint-oriented SV universe
and are merged with SURVIVOR. GRIDSS and CNVpytor are intentionally not added
to the SURVIVOR caller count:

- **Manta + DELLY**: paired-end, split-read and local-assembly evidence used for
  the primary SV event universe.
- **GRIDSS2**: independent breakpoint assembly. GRIDSS fundamentally reports
  breakpoint claims (BND representation), so its evidence is attached by
  breakpoint concordance rather than treated as another equivalent event caller.
- **CNVpytor**: independent read-depth CNV detection. Genome-wide DeepVariant
  SNPs are also imported to calculate BAF signals. CNVpytor DEL/DUP calls are
  matched to primary SVs by reciprocal overlap; unmatched calls are retained.
- **DeepVariant**: genome-wide WGS SNV/indel calling. Panel restriction is done
  downstream, not during discovery.
- **WhatsHap**: short-read local phasing/haplotagging. Interpret phase blocks as
  short-read evidence, not as equivalent to ONT long-range phasing.
- **MELTv2**: independent mobile-element insertion branch.
- **ExpansionHunter**: catalog-driven repeat-expansion genotyping.
- **Mutserve2**: mitochondrial homoplasmic/heteroplasmic SNV calling from chrM.
- **AnnotSV + VEP + Monarch**: annotation, panel/non-panel separation, phenotype
  context and candidate ranking, retaining the same conceptual goals as the LRS
  workflow.

## Why CNVpytor and GRIDSS are not merged into caller count

A CNVpytor duplication means a copy-number/read-depth change. A GRIDSS BND means
a breakpoint connection. Those are not equivalent claims. The final table
therefore contains explicit fields such as:

- `CNVPYTOR_RD_MATCH`
- `CNVPYTOR_RD_LEVEL`
- `CNVPYTOR_RECIPROCAL_OVERLAP`
- `GRIDSS_BREAKPOINT_MATCH`
- `GRIDSS_EVENT`
- `GRIDSS_MAX_QUAL`

The original Manta/DELLY agreement is labelled as a multi-caller-supported
companion callset rather than being interpreted as proof of high confidence.

## Main outputs

For each sample:

- `gene_discovery/<sample>_SRS_evidence_enhanced.tsv`
  - primary Manta/DELLY SV-gene rows plus independent GRIDSS and CNVpytor evidence.
- `gene_discovery/<sample>_CNVpytor_only_candidates.tsv`
  - read-depth CNVs not represented in the primary breakpoint SV universe.
- `gene_discovery/<sample>_GRIDSS_only_candidates.tsv`
  - GRIDSS breakpoint events not represented in the primary event universe.
- `mtdna/<sample>.mtDNA.vcf.gz`
  - Mutserve mitochondrial variants, including heteroplasmy according to the
    configured threshold.
- `mtdna/<sample>.mtDNA.vep.txt`
  - VEP consequence annotation of mitochondrial variants.
- `qc/<sample>.wgs_qc_summary.tsv`
  - combined samtools + mosdepth summary.
- `cnv/cnvpytor/<sample>.cnvpytor.pytor`
  - CNVpytor read-depth/SNP/BAF project file.

## Important setup checks

GRIDSS requires the FASTA to match the BAM reference and, with its default BWA
backend, requires the BWA index files next to the FASTA. The preflight rule
checks these when GRIDSS is enabled.

CNVpytor may require its human reference resources to be installed once in the
environment (the upstream command is `cnvpytor -download`) before the first
production run. Do this deliberately on a server with network access rather
than hiding a reference download inside Snakemake.

For GRIDSS, an hg38 ENCODE blacklist can be supplied in
`gridss_blacklist`. The workflow does not invent a local path; if no validated
blacklist file is installed, the value stays `null`.

`gridss_skip_softclip_realignment` is false by default. Enable it only when
the input alignments are known to contain appropriate split-read SA tags (for
example standard BWA-MEM output).

## Run

From this directory:

```bash
snakemake \
  -s Snakefile_SRS_WGS \
  --configfile config_srs_wgs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32
```

Before a full run, inspect the DAG:

```bash
snakemake \
  -s Snakefile_SRS_WGS \
  --configfile config_srs_wgs.yaml \
  --use-conda \
  --cores 32 \
  -n -p
```

## Methods grounding

The implementation follows the documented behavior of the upstream tools rather
than copying GATK-SV:

- GRIDSS2: Cameron et al., Genome Biology 2021; GRIDSS documentation.
- CNVpytor: Suvakov et al., GigaScience 2021.
- ExpansionHunter: Dolzhenko et al., Genome Research 2017; Bioinformatics 2019.
- Mutserve2 / mtDNA-Server methodology for mitochondrial heteroplasmy calling.
- DeepVariant for WGS SNV/indel calling.
- AnnotSV for structural-variant annotation.

The workflow intentionally keeps the evidence types separate so that technical
support can be interpreted biologically rather than reduced to a raw caller
count.


## SRS-specific thesis postprocessing

The separate `Snakefile_SRS_WGS_postprocess` mirrors only the LRS downstream
analyses whose assumptions remain valid for Illumina WGS:

- MitoCarta gene/pathway context;
- mechanism-aware SV-gene event ranking;
- optional exact gnomAD-SV site matching;
- genome-wide DeepVariant/VEP small-variant + SV same-gene review;
- short-read Samplot read-level review;
- cohort SV recurrence.

It intentionally does **not** copy LRS methylation, TLDR, Straglr, long-read
minimum-read-length Samplot settings, or long-range SV phasing assumptions.

Run after or together with the main SRS DAG:

```bash
snakemake \
  -s Snakefile_SRS_WGS_postprocess \
  --configfile config_srs_wgs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  all_srs_postprocess
```

The MitoCarta workbook and pathway GMX are deliberately configured as local
resources. They must exist at the configured paths. The repository currently
contains only the MitoCarta resource README, not the Broad data files
themselves. The workflow therefore fails visibly when the resources are absent
instead of silently substituting another dataset.

## Interpretation limits

`CALLER_COUNT` remains Manta/DELLY concordance only. A positive
`GRIDSS_BREAKPOINT_MATCH` is independent breakpoint-assembly evidence, and a
positive `CNVPYTOR_RD_MATCH` is independent read-depth evidence. These fields
must not be arithmetically added as if all algorithms represented the same
measurement.

CNVpytor BAF is calculated and stored in the `.pytor` project for regional
review, but the current enhanced TSV does not convert BAF into a binary
pathogenicity/support flag. That is intentional: BAF interpretation depends on
local SNP density, ploidy and event class and should not be reduced to a
generic threshold without validation.

The SRS SNV-SV table does not infer cis/trans between a small variant and an SV.
WhatsHap can phase short-read SNVs locally, but the breakpoint SV is not
long-range phased by this workflow.
