# SRS WGS pipeline — Illumina

This folder contains the **short-read whole-genome sequencing pipeline** used for
Illumina data.

The pipeline is kept separate from the LRS workflow. Nothing in the LRS folder
needs to be changed to run this SRS analysis.

The two files you normally need are:

- `Snakefile_SRS_WGS`
- `config_srs_wgs.yaml`

The `scripts/` folder contains only helper scripts used by this SRS workflow.
The `envs/` folder contains only the additional environments needed here.

---

## What the pipeline is trying to answer

The analysis is designed for hereditary optic neuropathy WGS and follows the
same general reasoning as the LRS pipeline:

1. Find structural variants across the genome.
2. Keep both known optic-neuropathy genes and non-panel genes.
3. Do not remove a candidate simply because only one type of evidence is
   available.
4. Add independent evidence when possible.
5. Keep mitochondrial biology separate from direct mtDNA variant calling.
6. Produce a final table that can be inspected biologically, not just a list of
   variants.

The optic-neuropathy panel is therefore used for **interpretation**, not to
restrict genome-wide discovery.

---

## Pipeline flow

```text
Illumina WGS BAM
      │
      ├── Quality control
      │     ├── mosdepth
      │     └── samtools flagstat / stats
      │
      ├── Structural variants from breakpoint evidence
      │     ├── Manta
      │     └── DELLY
      │             │
      │             └── SURVIVOR
      │                  → merged breakpoint SV callset
      │
      ├── Independent breakpoint / assembly support
      │     └── GRIDSS
      │
      ├── Genome-wide SNVs and indels
      │     └── DeepVariant
      │           │
      │           ├── WhatsHap local phasing
      │           └── SNP information for CNVpytor BAF
      │
      ├── Copy-number variants
      │     └── CNVpytor
      │           ├── read depth
      │           └── B-allele frequency context
      │
      ├── Combine SV + CNV evidence
      │     └── integrated SV VCF
      │
      ├── Annotation
      │     ├── AnnotSV
      │     └── VEP
      │
      ├── Gene interpretation
      │     ├── optic-neuropathy panel genes
      │     ├── non-panel genes
      │     ├── Monarch/HPO context
      │     └── optional MitoCarta context
      │
      ├── Other variant classes
      │     ├── ExpansionHunter → repeat expansions
      │     ├── MELT → mobile-element insertions
      │     └── Mutserve2 → mtDNA SNV / heteroplasmy
      │
      └── Diagnostic-support review
            ├── nuclear mitochondrial/bioenergetic candidates first
            ├── known disease genes
            ├── genome-wide candidates
            ├── secondary mtDNA branch
            ├── Samplot images/manifest
            └── optional candidate and Truvari benchmarks
```

---

## What each tool contributes

### Manta

Manta detects structural variants using paired-end, split-read and local
assembly evidence. It is one of the two main breakpoint-based SV callers.

### DELLY

DELLY is the second breakpoint-based caller. It provides paired-end and
split-read evidence that can be compared with Manta.

### SURVIVOR

SURVIVOR merges the Manta and DELLY calls into one breakpoint-based SV set.

The file containing variants detected by at least two callers is named
`multi_caller`, not `high_confidence`.

This is deliberate: two callers can use related read evidence, so agreement
between Manta and DELLY is useful support but does not by itself prove that an
SV is real.

### GRIDSS

GRIDSS is kept separate from SURVIVOR.

GRIDSS performs breakpoint assembly and represents many rearrangements as
breakends. The pipeline therefore uses GRIDSS as **independent supporting
evidence**.

For events with two defined breakpoints, the final table distinguishes:

- `YES_BOTH_BREAKPOINTS`
- `PARTIAL_ONE_BREAKPOINT`
- `NO`

A one-breakpoint GRIDSS match is not treated as equivalent to support for the
complete rearrangement.

### DeepVariant

DeepVariant is run genome-wide.

This is important because the project is WGS-based. Small-variant discovery is
therefore not limited to the known optic-neuropathy panel.

The DeepVariant VCF is also used for:

- local WhatsHap phasing;
- SNP/allelic information used by CNVpytor.

The complete genome-wide VCF is preserved. A second PASS-focused VCF and VEP
table are generated for loci in `nuclear_mito_candidate_bed`, so SNV/indel
candidates in the main nuclear mitochondrial/bioenergetic hypothesis are easy
to review without pretending that the rest of the WGS callset was never made.
The repository default is the curated optic-neuropathy mitochondrial-gene BED;
replace it with a broader validated nuclear-mitochondrial BED when appropriate.

### CNVpytor

CNVpytor adds information that the breakpoint callers do not provide directly:
**read-depth copy-number evidence**.

The workflow follows the CNVpytor documented sequence:

```text
BAM
 ↓
read depth
 ↓
histograms
 ↓
partitioning
 ↓
CNV calls

DeepVariant VCF
 ↓
SNP import
 ↓
SNP mask
 ↓
BAF calculation
```

Two bin sizes are retained by default:

- 10 kb: main CNV analysis;
- 100 kb: large-CNV context.

The file `*.reference_info.txt` records the CNVpytor `-ls` output so you can
check that hg38 and the required GC/mask resources were recognized.

The pipeline does **not** use CNVpytor's prototype combined RD+BAF caller as the
main CNV callset. Standard read-depth calls are used for discovery; BAF remains
supporting information.

### Combining SV and CNV evidence

The script:

`scripts/combine_sv_cnv_evidence.py`

does three simple things:

1. keeps the Manta+DELLY breakpoint events;
2. adds CNVpytor read-depth evidence to overlapping DEL/DUP events;
3. keeps CNVpytor-only DEL/DUP calls if they pass the configured CNVpytor
   quality checks.

The resulting file is:

```text
<sample>/sv/integrated/<sample>_integrated_SV.vcf.gz
```

This is the structural-variant file that continues to AnnotSV and VEP.

### VEP

VEP is run with `--flag_pick`, not `--pick`.

Therefore all transcript consequences are retained. The script:

`scripts/summarize_vep_transcripts.py`

creates a readable per-SV/per-gene summary without deleting the other
transcripts.

### ExpansionHunter

ExpansionHunter is used for the repeat-expansion catalog.

It is a separate variant class and is **not** counted as another generic SV
caller.

### MELT

MELT is the optional mobile-element insertion branch.

It remains disabled until the local MELTv2 installation and reference files are
configured.

### Mutserve2

Mutserve2 is used for mtDNA SNVs and heteroplasmy.

This is included because hereditary optic neuropathies can be caused directly by
mitochondrial-genome variants.

The output is summarized by:

`scripts/summarize_mtdna_variants.py`

The current branch is deliberately described as **mtDNA SNV/heteroplasmy
analysis**, not complete mtDNA variant detection, because Mutserve2 VCF output
does not provide a complete mtDNA indel solution.

### MitoCarta

MitoCarta answers a different question from Mutserve2.

Mutserve2 evaluates variants **inside the mitochondrial genome**.

MitoCarta is used to ask whether a **nuclear gene affected by an SV encodes a
mitochondrial protein or belongs to a mitochondrial pathway**.

This is particularly useful for the optic-neuropathy hypothesis.

---

## Final technical-support labels

The final SRS table uses readable labels:

```text
MULTI_CALLER_PLUS_READ_DEPTH
MULTI_CALLER_PLUS_GRIDSS
MULTI_CALLER
SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE
SINGLE_CALLER
READ_DEPTH_ONLY_CNV
REVIEW_REQUIRED
```

These labels describe only the technical evidence supporting the call.

They are **not** pathogenicity classifications and they are **not** ACMG
classes.

---

## Files you should edit first

Open:

```text
config_srs_wgs.yaml
```

and change the real server paths for:

```yaml
samples:
path:
ref:
candidate_genes_list:
nuclear_mito_candidate_bed:
exclude_bed:
expansionhunter_catalog:
vep_cache_dir:
annotsv_annotations_dir:
monarch_nodes:
monarch_edges:
```

Optional tools/resources such as MELT and MitoCarta can remain disabled until
their local files are ready.

For the intended mitochondrial-biology analysis, install the MitoCarta 3.0
inventory and pathway GMX, set their paths, and change `mitocarta_enabled` to
`true`. This promotes nuclear-encoded mitochondrial genes and preserves their
pathway/subcompartment context. Mutserve remains a separate secondary mtDNA
branch.

## Install and validate the environments

Create the small runner environment once:

```bash
mamba env create -f environment.yaml
mamba activate srs-wgs-runner
```

After editing `config_srs_wgs.yaml`, resolve and create every environment used
by the selected DAG without running an analysis:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --conda-create-envs-only --cores 1
```

Then validate paths, indexes, optional-resource switches, and the DAG:

```bash
snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --cores 1 reference_checks/srs_preflight.tsv

snakemake -s Snakefile_SRS_WGS --configfile config_srs_wgs.yaml \
  --use-conda --cores 32 --dry-run
```

Preflight stops before expensive callers if a BAM/index, FASTA index, GRIDSS
BWA index, annotation bundle, HPO/Monarch resource, or enabled optional resource
is missing. Each sample also receives `qc/*.wgs_qc_summary.tsv`; threshold
failures are marked `REVIEW` instead of silently discarding calls.

## Candidate review and benchmarking

The principal review output is:

```text
<sample>/diagnostic_review/<sample>.diagnostic_candidates.tsv
```

It is an ordered research/diagnostic-support table, not a pathogenicity
classification. Samplot writes a manifest beside it. Interchromosomal BNDs are
explicitly retained for manual two-breakpoint review instead of being drawn as
false same-chromosome intervals.

To benchmark prioritization, copy `validation/candidate_truth.template.tsv`,
enter known positive-control genes/IDs, and set `candidate_truth_tsv`. To
benchmark SV detection against a sample truth VCF, add that sample under
`sv_truth_vcfs`; Truvari output is written under `<sample>/benchmark/truari/`.

---

## Dry run

From this folder:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main_v3/snakemake_pipelines/srs_wgs_pipeline
```

first run:

```bash
snakemake \
    -s Snakefile_SRS_WGS \
    --configfile config_srs_wgs.yaml \
    --use-conda \
    --conda-prefix "/home/casadei7/snakemake_envs/envs/" \
    --cores 32 \
    -n
```

Only after the dry run is clean should the complete workflow be launched.

---

## Main methodological references

The choices in this workflow follow the published/documented methods of the
tools rather than treating them as interchangeable callers:

- Manta — Chen et al., *Bioinformatics*, 2016.
- DELLY — Rausch et al., *Bioinformatics*, 2012.
- GRIDSS2 — Cameron et al., *Genome Biology*, 2021.
- CNVpytor — Suvakov et al., *GigaScience*, 2021.
- DeepVariant — Poplin et al., *Nature Biotechnology*, 2018.
- ExpansionHunter — Dolzhenko et al., *Genome Research*, 2017 and
  *Bioinformatics*, 2019.
- Mutserve / mtDNA-Server — Weissensteiner et al. and subsequent mtDNA-Server
  developments.

---

## Important interpretation limits

The workflow intentionally does not make conclusions that the available data
cannot support.

- Multi-caller support is not equivalent to variant truth.
- GRIDSS partial breakpoint support is kept separate from full breakpoint
  support.
- CNVpytor read depth supports copy-number change, not the exact breakpoint.
- Short-read phasing is not treated as equivalent to long-read phasing.
- ExpansionHunter and MELT are separate variant classes, not extra votes for all
  SVs.
- mtDNA SNV analysis and nuclear mitochondrial-gene annotation are kept
  separate.
- Family inheritance is not inferred when trio/family data are unavailable.
