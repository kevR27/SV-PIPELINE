# SV-PIPELINE

Long- and short-read WGS workflows for genome-wide structural-variant discovery with an optic-neuropathy/rare-neurological-disease interpretation layer.

## Active LRS workflow

Use:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main/snakemake_pipelines/lrs

snakemake \
  --snakefile Snakefile_LRS_update \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --printshellcmds
```

`Snakefile_LRS_update` reads `config_lrs.yaml` from the same directory.

Before a full run, inspect the DAG/dry run:

```bash
snakemake \
  --snakefile Snakefile_LRS_update \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --dry-run \
  --printshellcmds
```

## LRS biological/data-flow design

```text
ONT BAM
│
├── mosdepth                           sequencing-depth QC
├── Clair3 -> WhatsHap                complementary panel SNV/indel phasing
│
├── Sniffles2 ─┐
├── cuteSV     ├─ caller evidence filter -> Jasmine -> MASTER SV VCF
└── DELLY      ┘                              │
                                              ├── caller-support summary
                                              ├── AnnotSV genome-wide
                                              │    └── panel-only derived view
                                              ├── VEP supplementary annotation
                                              ├── Monarch/HPO gene discovery
                                              └── integrated SV/gene evidence TSV

BAM -> dedicated Sniffles2 v2.6.2 query -> needLR
                                           └── population-AF evidence joined
                                               back to the master SV table

BAM -> Straglr                         tandem-repeat expansion discovery
BAM -> TLDR                            mobile-element insertion confirmation
BAM -> modkit                          methylation calls
Clair3 + filtered Sniffles2 -> LongPhase  same-patient SNP/SV phasing
```

### Important interpretation rules

- The complete sorted/indexed Jasmine VCF is the **master per-patient SV universe**.
- The `>=2 callers` VCF is a **companion evidence set**, not a replacement for single-caller variants.
- AnnotSV and VEP operate on the complete Jasmine master VCF.
- needLR is used as a **supplementary ONT population-frequency annotation branch**. It is not allowed to remove variants from the master callset.
- BNDs and SVs outside needLR's supported size range remain in the integrated table and are labelled as not evaluable by needLR.
- SV calling is genome-wide; the optic-neuropathy panel is applied as a downstream annotation/view.
- The Monarch/HPO ranking is a discovery-priority analysis, not an automatic pathogenicity classification.

## needLR compatibility

The dedicated needLR query is generated with Sniffles2 2.6.2, matching the caller version used for the needLR v4 ONT reference as closely as possible. The discovery Sniffles2 call remains separate so clinically useful discovery settings do not have to be forced onto the population-reference comparison.

## Legacy helper scripts

`scripts/prepare_needlr_vcf.py` and `scripts/fix_needlrvcf_for_annotsv.py` are retained for manual troubleshooting/backward compatibility, but they are intentionally **not part of the active default DAG**. AnnotSV no longer depends on a needLR-derived VCF.

## Active SRS workflow

Use:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main/snakemake_pipelines/srs

snakemake \
  --snakefile Snakefile_SRS \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --printshellcmds
```

Always inspect the DAG before a full run:

```bash
snakemake \
  --snakefile Snakefile_SRS \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --dry-run \
  --printshellcmds
```

### SRS biological/data-flow design

```text
Illumina BAM
│
├── mosdepth                              sequencing-depth QC
├── DeepVariant -> WhatsHap              complementary panel SNV/indel arm
│
├── Manta ─┐
└── DELLY ─┴─ caller evidence filter -> SURVIVOR -> MASTER SV VCF
                                                  │
                                                  ├── caller-support summary
                                                  ├── AnnotSV genome-wide
                                                  │    └── panel-only derived view
                                                  ├── VEP supplementary annotation
                                                  ├── Monarch/HPO gene discovery
                                                  └── integrated SV/gene evidence TSV

BAM -> ExpansionHunter                   catalog-driven STR genotyping
BAM -> MELTv2 Single                     independent MEI discovery (optional)
```

### SRS interpretation rules

- The complete sorted/indexed SURVIVOR VCF is the **master per-patient SRS SV universe**.
- The `>=2 callers` VCF is a **companion evidence set**, not a replacement for single-caller calls.
- Manta and DELLY call genome-wide. The optic-neuropathy panel is applied only to downstream annotation/prioritization; DeepVariant remains panel-targeted by design.
- Manta support is normalized from alternate paired-read plus split-read counts (`PR_alt + SR_alt`). DELLY uses `PE + SR`, with `DV + RV` as a fallback.
- AnnotSV and VEP use the complete SURVIVOR master VCF. The panel AnnotSV file is derived from the full annotation rather than rerunning a restricted annotation.
- needLR is not used for SRS. The shared integrated-table schema reports `NEEDLR_STATUS=NOT_APPLICABLE_SRS`.
- ExpansionHunter and MELT remain independent evidence branches; they are not silently merged into the breakpoint-SV master VCF. A companion `*_integrated_SV_gene_with_orthogonal_evidence.tsv` attaches coordinate-aware MELT matches and ExpansionHunter locus overlaps without changing caller counts or removing rows.
- MELT is disabled until its manually licensed `MELT.jar`, gene BED, and family ZIP paths are configured. Run one transposon family per isolated working directory.
- Candidate rankings are discovery priorities, not automatic pathogenicity classifications.

Generate the applicable per-sample thesis plots with the SRS caller order:

```bash
python ../../plots/run_thesis_plots.py \
  --root /DATA/srs_results \
  --sample patient01 \
  --platform srs
```
