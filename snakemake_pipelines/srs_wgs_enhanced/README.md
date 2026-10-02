# Enhanced SRS WGS pipeline (Illumina)

This folder is an **isolated SRS-only workflow**. It was added specifically so that the
existing LRS Snakefiles, LRS post-processing workflow, and shared LRS-sensitive scripts
do not need to be edited.

Main entry points:

- `Snakefile_SRS_enhanced`
- `config_srs_enhanced.yaml`
- `scripts/` — helper scripts used only by this enhanced SRS workflow
- `envs/` — environments used only by this enhanced SRS workflow

The workflow still reads the existing repository panel/reference files and shared
annotation scripts, but it does not modify them.

## Biological design

The workflow deliberately separates **variant discovery modalities** instead of
treating every tool as an interchangeable SV caller.

### 1. Breakpoint discovery: Manta + DELLY

Manta and DELLY remain the primary symbolic breakpoint callers. Their filtered calls
are merged with SURVIVOR to create the breakpoint-event universe.

`SUPP >= 2` is now labelled **caller-concordant**, not “high confidence”. Manta and
DELLY can be driven by related paired-end/split-read evidence, so caller agreement is
supporting evidence rather than an independent proof of validity.

### 2. GRIDSS: independent assembly/breakend evidence

GRIDSS is run as a separate assembly/breakend branch. It is **not inserted as a third
SURVIVOR caller**, because GRIDSS represents rearrangements primarily as BND records
and performs local assembly. Its breakpoint evidence is attached downstream to master
events by proximity.

GRIDSS documentation recommends an ENCODE exclusion list for human data. For GRCh38
with UCSC chromosome naming, configure the local ENCFF356LFX BED in
`gridss_exclude_bed` when available.

The GRIDSS documentation also reports assembly memory requirements around 31 GB for
its standard parallel assembly stage; the rule therefore requests 32 GB and uses up to
8 threads by default.

### 3. DeepVariant: whole-genome SNV/indel calling

DeepVariant is now run **genome-wide**, not on the optic-neuropathy BED. Panel and
non-panel interpretation happens after calling. This keeps WGS discovery symmetric:
genome-wide SV discovery is no longer paired with panel-only small-variant calling.

DeepVariant output also provides the SNP information used by CNVpytor for BAF context
and by WhatsHap for local short-read phasing.

### 4. CNVpytor: read-depth CNV + BAF context

CNVpytor adds an independent read-depth modality that Manta/DELLY do not replace.

The implemented workflow follows the documented sequence:

1. import read depth from BAM;
2. calculate histograms;
3. partition read depth;
4. import SNPs from the DeepVariant VCF;
5. apply the SNP mask;
6. calculate BAF histograms;
7. call CNVs at two resolutions.

The primary default is 10 kb and the secondary context scale is 100 kb. The workflow
does **not** rely on CNVpytor's combined RD+BAF caller as the master CNV caller because
the project documentation labels that combined caller as prototype. RD CNV calls are
the discovery signal; BAF remains an orthogonal context/evidence layer.

For unmatched CNVpytor calls to enter the enhanced master universe, the default QC
uses documented CNVpytor fields:

- q0 <= 0.5
- pN <= 0.5
- e-val1 <= 1e-4

These are research filters, not pathogenicity criteria.

Same-type DEL/DUP events are attached to Manta/DELLY calls when minimum reciprocal
overlap is >= 0.5. A passing CNVpytor call without a breakpoint match is retained as a
`CNVPYTOR_ONLY` symbolic DEL/DUP so a depth-only CNV is not silently discarded.

### 5. Enhanced master SV universe

The master VCF contains:

- Manta/DELLY breakpoint-union events;
- CNVpytor RD evidence attached to matching DEL/DUP events;
- passing CNVpytor-only DEL/DUP events.

GRIDSS remains orthogonal evidence rather than redefining these events.

This makes the evidence model explicit:

- breakpoint evidence;
- independent read-depth evidence;
- BAF context;
- assembly/breakend evidence;
- specialized MEI and STR evidence.

### 6. ExpansionHunter

ExpansionHunter remains a catalog-driven repeat-expansion genotyper. It is a separate
variant class and is **not counted as generic SV-caller support**.

Its own documentation states that it searches BAM/CRAM reads spanning, flanking, or
contained within catalogued repeats.

### 7. MELTv2

MELT remains an optional independent mobile-element-insertion branch. It is disabled
by default because a valid MELTv2 installation and transposon resources must be
configured locally.

### 8. mtDNA: Mutserve2

Because hereditary optic neuropathy can be caused by mitochondrial-genome variants,
the enhanced workflow adds an mtDNA branch.

The pipeline extracts primary non-duplicate chrM alignments and runs Mutserve2 for
homoplasmic/heteroplasmic mtDNA variants. The mtDNA reference is extracted directly
from the same hg38 FASTA used for the BAM, avoiding contig/reference mismatch.

Important limitation: Mutserve2 documentation states that its VCF output currently
does not include indels. Therefore this branch is explicitly treated as
**mtDNA SNV/heteroplasmy evidence**, not a complete mtDNA indel workflow. The pipeline
does not hide this limitation.

### 9. Annotation and gene discovery

The enhanced master VCF is passed to the same conceptual annotation framework used in
the project:

- AnnotSV genome-wide annotation;
- panel-only derivative table;
- VEP supplementary annotation;
- panel/non-panel gene separation;
- phenotype/HPO context;
- allele-level research assessment.

Optional MitoCarta context can be enabled once the local MitoCarta3.0 workbook and
MitoPathways GMX files are configured.

### 10. Technical evidence classes

The final SRS augmentation produces descriptive classes such as:

- `BREAKPOINT_CONCORDANT_PLUS_RD`
- `BREAKPOINT_CONCORDANT_PLUS_ASSEMBLY`
- `BREAKPOINT_CALLER_CONCORDANT`
- `SINGLE_BREAKPOINT_CALLER_PLUS_ORTHOGONAL`
- `SINGLE_BREAKPOINT_CALLER`
- `RD_ONLY_CNV`
- `REVIEW_REQUIRED`

These labels describe the **technical evidence architecture only**. They are not ACMG
classification and must not be interpreted as pathogenic/benign labels.

## Why GATK-SV is not embedded

This pipeline borrows the *principle* of combining different evidence modalities, but
does not embed GATK-SV. The goal is a transparent Snakemake workflow suitable for the
existing server and thesis analysis while preserving the current project structure.

The design instead combines tools that can run directly on a local Linux server:

- Manta / DELLY: breakpoint calling
- GRIDSS: assembly/breakend evidence
- CNVpytor: read depth + BAF context
- DeepVariant: small variants
- ExpansionHunter: repeat expansions
- MELT: MEIs
- Mutserve2: mtDNA SNV/heteroplasmy

## Run

Edit `config_srs_enhanced.yaml` first, especially BAM paths, result path, reference,
AnnotSV, VEP, Monarch resources, and optional GRIDSS/MELT/MitoCarta resources.

Then:

```bash
snakemake \
  -s Snakefile_SRS_enhanced \
  --configfile config_srs_enhanced.yaml \
  --use-conda \
  --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32 \
  -n
```

Run the dry-run first. If the DAG is correct:

```bash
snakemake \
  -s Snakefile_SRS_enhanced \
  --configfile config_srs_enhanced.yaml \
  --use-conda \
  --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
  --cores 32
```

## Methods used as design references

- Manta: Chen et al., *Bioinformatics* 2016; paired-end/split-read/local assembly SV
  discovery.
- DELLY: Rausch et al., *Bioinformatics* 2012 and subsequent DELLY2 development.
- GRIDSS2: Cameron et al., *Genome Biology* 2021; genome-wide breakend assembly,
  split-read and read-pair evidence.
- CNVpytor: Suvakov et al., *GigaScience* 2021; CNV/CNA analysis from read depth and
  allele imbalance.
- DeepVariant: Poplin et al., *Nature Biotechnology* 2018 and current Google
  DeepVariant WGS implementation.
- ExpansionHunter: Dolzhenko et al., *Genome Research* 2017; *Bioinformatics* 2019.
- Mutserve / mtDNA-Server: Weissensteiner et al., *NAR* 2016; *Genome Research* 2021;
  mtDNA-Server 2, *NAR* 2024.

## Deliberate exclusions / cautions

- The CNVpytor combined RD+BAF caller is not used as the master caller because its
  documentation currently labels it prototype.
- GRIDSS is not merged blindly with symbolic Manta/DELLY records.
- Caller count is not equated with variant truth.
- Short-read WhatsHap output is local phasing evidence; it is not treated as
  equivalent to long-read phasing.
- Mutserve VCF is not described as a complete mtDNA indel solution.
- No trio/inheritance inference is invented when family data are absent.
