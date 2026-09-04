# SV-PIPELINE

A Snakemake-based structural-variant analysis pipeline for patients with undiagnosed rare neurological disease and suspected optic neuropathy.

The repository provides two complementary analysis workflows:

- **SRS** — Short-Read Sequencing
- **LRS** — Long-Read Sequencing, primarily Oxford Nanopore

Both workflows preserve sequencing-platform-specific variant discovery while converging on a common downstream clinical interpretation strategy based on **AnnotSV, gene extraction, human gene–phenotype associations, and candidate ranking**.

---

## 1. Concept

The main design principle is:

```text
                         SV-PIPELINE
                              │
               ┌──────────────┴──────────────┐
               │                             │
              SRS                           LRS
        Short-read BAM                 Long-read BAM
               │                             │
        Manta + Delly                Sniffles2 + cuteSV
               │                             │
           SURVIVOR                       Jasmine
               │                             │
               │                          needLR
               │                             │
               └──────────────┬──────────────┘
                              │
                           AnnotSV
                              │
                 ┌────────────┴────────────┐
                 │                         │
          Panel-focused                Genome-wide
            analysis                    annotation
                 │                         │
                 │                  SV-associated genes
                 │                         │
                 │                  Remove panel genes
                 │                         │
                 │                  Human HGNC/HPO layer
                 │                         │
                 │                  Candidate ranking
                 │
        Platform-specific analyses
```

The **genome-wide AnnotSV result give us an overview of all genes that overlaps from structural variants discovery obtained through SV callers**; while the panel-only AnnotSV result is used as a separate clinical view that uses all known genes that are directly associated to Optic Neuropathies 

---

# 2. Repository structure

```text
SV-PIPELINE/
│
├── PANEL_OA/
│   ├── optic_neuropathy_panel.bed
│   └── optic_neuropathy_genes.txt
│
├── envs/
│   ├── annotsv.yaml
│   ├── clair3.yaml
│   ├── cutesv.yaml
│   ├── deepvariant.yaml
│   ├── delly.yaml
│   ├── expansionhunter.yaml
│   ├── jasmine.yaml
│   ├── longphase.yaml
│   ├── manta.yaml
│   ├── melt.yaml
│   ├── modkit.yaml
│   ├── monarch.yaml
│   ├── mosdepth.yaml
│   ├── needlr.yaml
│   ├── repeatmasker.yaml
│   ├── sniffles2.yaml
│   ├── straglr.yaml
│   ├── stranger.yaml
│   ├── survivor.yaml
│   ├── tldr.yaml
│   ├── vep.yaml
│   └── whatshap.yaml
│
├── reference/
│   ├── optic_neuropathy_hpo_terms.tsv
│   ├── variant_catalog_hg38_repeats.json
│   └── [additional local/reference resources]
│
├── scripts/
│   ├── extract_annotsv_genes.py
│   ├── extract_nonpanel_genes.py
│   ├── monarch_human_gene_phenotypes.py
│   └── rank_sv_gene_candidates.py
│
├── snakemake_pipelines/
│   ├── srs/
│   │   ├── Snakefile_SRS
│   │   └── config_srs.yaml
│   │
│   └── lrs/
│       ├── Snakefile_LRS
│       └── config_lrs.yaml
│
├── README.md
└── .gitattributes
```
---

# 3. SRS pipeline

## 3.1 Purpose

The SRS workflow is designed for BAM files obtained from short-read sequencing such as Illumina which uses pair end sequencing mode to obtained the whole genome sequence of patients 

### Main workflow

```text
Sorted BAM
   │
   ├── mosdepth
   │       └── genome coverage
   │
   ├── Manta
   │       └── SV calls
   │
   └── Delly
           └── SV calls
               │
               ▼
            SURVIVOR
               │
               ▼
       Merged genome-wide SV VCF
               │
        ┌──────┴─────────┐
        │                │
     AnnotSV          VEP
        │                │
        │         full + panel view
        │
        ├── panel-only AnnotSV
        │
        └── gene discovery
                │
        ┌───────┴────────┐
        │                │
     SV genes       non-panel genes
                         │
                  Monarch HGNC/HPO
                         │
                   ranked candidates
```

The current SRS Snakefile implements mosdepth for QC , Manta, + Delly for SV discovery , SURVIVOR, AnnotSV, VEP, MELT, ExpansionHunter, DeepVariant and WhatsHap.

## 3.2 SRS tools

| Tool | Function |
|---|---|
| **mosdepth** | Genome-level sequencing coverage/QC |
| **Manta** | Short-read structural variant discovery |
| **Delly** | Short-read structural variant discovery |
| **SURVIVOR** | Merges Manta and Delly SV callsets |
| **AnnotSV** | SV annotation and clinical/gene annotation |
| **VEP** | Variant Effect Predictor annotation |
| **MELT** | Mobile element insertion detection |
| **ExpansionHunter** | Short tandem-repeat expansion analysis |
| **DeepVariant** | Small-variant calling |
| **WhatsHap** | Variant phasing and haplotagging |

### SV merging

Manta and Delly are intentionally retained as independent callers. SURVIVOR creates the combined SRS SV callset:

```text
Manta ──┐
        ├──> SURVIVOR ──> merged_SV.vcf
Delly ──┘
```

A separate **SUPP>=2 high-confidence companion SV set** is also defined for review, without replacing the complete merged callset.

---

# 4. LRS pipeline

## 4.1 Purpose

The LRS workflow is designed primarily for long-read sequencing data, particularly Oxford Nanopore BAM files.
[INSERT PIPELINE OF BASECALLING FROM POD5 USING DORADO]

### Main workflow

```text
Long-read BAM
    │
    ├── mosdepth
    │
    ├── Clair3
    │      └── small variants
    │
    ├── Sniffles2
    │
    └── cuteSV
           │
           ▼
        Jasmine
           │
           ▼
       merged SV VCF
           │
         needLR
           │
           ▼
        AnnotSV
           │
      ┌────┴───────────────┐
      │                    │
 Panel-only          Genome-wide
 AnnotSV              AnnotSV
                           │
                           ▼
                    SV-associated genes
                           │
                    remove panel genes
                           │
                    human HGNC/HPO
                           │
                    ranked candidates
```

The current LRS workflow preserves its platform-specific upstream analysis and then converges on the same downstream clinical/gene-discovery architecture used by SRS.

## 4.2 LRS tools

| Tool | Function |
|---|---|
| **mosdepth** | Coverage/QC |
| **Clair3** | Small-variant calling from ONT reads |
| **WhatsHap** | SNP phasing/haplotagging |
| **Sniffles2** | Long-read SV discovery |
| **cuteSV** | Long-read SV discovery |
| **Jasmine** | Merges long-read SV callsets |
| **needLR** | ONT population/reference SV annotation/filtering |
| **AnnotSV** | SV annotation |
| **VEP** | Transcript-level SV annotation |
| **RepeatMasker** | Transposable-element/insertion classification |
| **Straglr** | Repeat-expansion analysis |
| **modkit** | ONT methylation extraction |
| **LongPhase** | Joint SNP + SV phasing |
| **TLDR** | Confirmatory mobile-element insertion analysis |

The current LRS implementation uses Sniffles2 and cuteSV followed by Jasmine, then needLR and AnnotSV.

---

# 5. LRS-specific analyses

## RepeatMasker

Jasmine insertion calls are converted into insertion sequences and analyzed with RepeatMasker using:

```text
RepeatMasker -species human
```

The purpose is to characterize insertion sequences and support transposable-element/mechanism classification.

## Straglr

Straglr is available for repeat-expansion analysis using:

```text
BAM + reference + optic-neuropathy panel BED
```

It is configured through `straglr_script` rather than a hard-coded `/path/to/...` executable.

The Straglr rule remains available, but it is not part of the current default `rule all` target.

## needLR

needLR is used after Jasmine:

```text
Sniffles2
     +
cuteSV
     │
  Jasmine
     │
  needLR
```

The workflow collects the generated `*_RESULTS.tsv` and `*_RESULTS.vcf.gz` files and indexes the final VCF.

needLR should be understood here as a **population/reference SV annotation and filtering layer**, not as a primary SV discovery caller.

## modkit

The LRS workflow includes ONT methylation extraction with modkit.

The CpG region is configurable through:

```yaml
cpg:
```

The current default configuration sets:

```yaml
cpg: null
```

thus, no chromosome region is imposed.

## LongPhase

LongPhase combines the Clair3 SNP callset with the Jasmine SV callset:

```text
Clair3 SNPs
     +
Jasmine SVs
     +
ONT BAM
     │
     ▼
 LongPhase
```

The LRS workflow explicitly invokes LongPhase with ONT mode enabled.

## TLDR

TLDR provides confirmatory mobile-element insertion analysis using:

```text
ONT BAM + reference + teref.ont.human.fa
```

and writes a per-sample TLDR table.

---

# 6. AnnotSV architecture

AnnotSV is deliberately run in two modes.

## Genome-wide annotation

```text
SV VCF
  │
AnnotSV
  │
candidateGenesFiltering 0
  │
complete genome-wide annotation
```

The complete result is used as the source for novel/non-panel gene discovery.

## Panel-only annotation

```text
SV VCF
  │
AnnotSV
  │
candidateGenesFiltering 1
  │
optic-neuropathy panel
```

The panel-only result is a focused clinical view and does not replace the genome-wide annotation.

Both workflows also support the configurable AnnotSV benign allele-frequency threshold:

```yaml
annotsv_benign_af: 0.01
```

The current workflows pass this value to AnnotSV as `-benignAF`.

---

# 7. Genome-wide gene discovery

This is a key component of the pipeline.

The objective is to identify genes associated with SVs **outside the predefined optic-neuropathy panel** and then provide human phenotype evidence that can help prioritize candidate genes.

## Step 1 — Extract SV-associated genes

```text
Genome-wide AnnotSV TSV
          │
          ▼
extract_annotsv_genes.py
          │
          ▼
SV_genes.txt
```

The script extracts unique gene symbols from AnnotSV's `Gene_name` field.

## Step 2 — Remove genes already present in the panel

```text
SV_genes.txt
      +
optic_neuropathy_genes.txt
      │
      ▼
extract_nonpanel_genes.py
      │
      ▼
nonpanel_genes.txt
```

This prevents the discovery stage from simply returning genes that were already included in the clinical panel.

## Step 3 — Human gene → phenotype associations

```text
nonpanel_genes.txt
       +
optic-neuropathy HPO anchors
       │
       ▼
Monarch Knowledge Graph
       │
       ▼
human_gene_phenotypes.tsv
```

The Monarch script deliberately resolves human genes through **HGNC entities** and retrieves human gene-to-phenotype associations. It does not traverse model-organism ortholog associations.

The included HPO terms are **reference optic-neuropathy phenotype anchors**.

Important:

> The HPO layer represents gene-associated phenotype evidence. It is not a prediction of the individual patient's actual phenotype,but from here, it will be possible to infer possible clinical diagnosis 

## Step 4 — Candidate ranking

```text
Genome-wide AnnotSV
        +
non-panel genes
        +
human gene/HPO associations
        +
optic-neuropathy panel
        │
        ▼
rank_sv_gene_candidates.py
        │
        ▼
ranked_candidates.tsv
```

This produces a prioritization table for genome-wide candidate genes.

The resulting ranking is **candidate prioritization, not an automatic pathogenicity classification**.

---

# 8. HPO strategy

The pipeline does not require a complete patient-specific HPO profile for genome-wide discovery.

The reference file is:

```text
reference/optic_neuropathy_hpo_terms.tsv
```

The current workflow uses reference phenotype anchors including:

```text
HP:0000648   Optic atrophy
HP:0007663   Reduced visual acuity
HP:0000543   Optic disc pallor
HP:0000649   Abnormality of visual evoked potentials
```

These terms are used to identify phenotype associations relevant to optic neuropathy at the **gene level**, rather than claiming that they describe the patient's observed phenotype.

---

# 9. Input requirements

Each workflow requires a sorted BAM and the appropriate configuration file.

## SRS

Edit:

```text
snakemake_pipelines/srs/config_srs.yaml
```

Example:

```yaml
samples:
  patient01: /path/to/patient01.sorted.bam
  patient02: /path/to/patient02.sorted.bam

path: "/path/to/srs_results/"
ref: "/path/to/hg38.fa"
threads: 16

bed: "PANEL_OA/optic_neuropathy_panel.bed"
candidate_genes_list: "PANEL_OA/optic_neuropathy_genes.txt"
```

Additional SRS resources include:

```yaml
exclude_bed:
expansionhunter_catalog:
melt_dir:
melt_transposons:
melt_genes_bed:
vep_cache_dir:
```

The values shown in the repository's current configuration are **examples/placeholders and must be adapted to the real analysis machine**.

## LRS

Edit:

```text
snakemake_pipelines/lrs/config_lrs.yaml
```

Example:

```yaml
samples:
  patient01: /path/to/patient01.sort.bam
  patient02: /path/to/patient02.sort.bam

path: "/path/to/lrs_results/"
ref: "/path/to/hg38.fa"
threads: 16

bed: "PANEL_OA/optic_neuropathy_panel.bed"
candidate_genes_list: "PANEL_OA/optic_neuropathy_genes.txt"

cpg: null

condaenv: "envs/"
```

Additional LRS machine-specific resources include:

```yaml
modkitenv:
clair3env:
vntr_bed:
vep_cache_dir:
needlr_backend_files:
tldr_elements:
```

These paths must point to actual resources on the analysis system.

---

# 10. Conda environments

Each major tool has its own environment definition.

Current environments:

```text
annotsv.yaml
clair3.yaml
cutesv.yaml
deepvariant.yaml
delly.yaml
expansionhunter.yaml
jasmine.yaml
longphase.yaml
manta.yaml
melt.yaml
modkit.yaml
monarch.yaml
mosdepth.yaml
needlr.yaml
repeatmasker.yaml
sniffles2.yaml
straglr.yaml
stranger.yaml
survivor.yaml
tldr.yaml
vep.yaml
whatshap.yaml
```

Snakemake selects the environment associated with each rule through the `conda:` directive. The repository currently provides environment definitions for the platform-specific and shared downstream analyses.

---

# 11. Running the SRS pipeline

From:

```text
snakemake_pipelines/srs/
```

activate the Snakemake environment:

```bash
conda activate snakemake
```

Then:

```bash
snakemake -s Snakefile_SRS --lint
```

Dry-run:

```bash
snakemake -s Snakefile_SRS -n --cores 1
```

Production execution:

```bash
snakemake -s Snakefile_SRS --cores 16 --use-conda
```

The working directory should be:

```text
SV-PIPELINE/snakemake_pipelines/srs
```

because the Snakefile refers to:

```python
configfile: "config_srs.yaml"
```

---

# 12. Running the LRS pipeline

From:

```text
snakemake_pipelines/lrs/
```

activate Snakemake:

```bash
conda activate snakemake
```

Lint:

```bash
snakemake -s Snakefile_LRS --lint
```

Dry-run:

```bash
snakemake -s Snakefile_LRS -n --cores 1
```

Production execution:

```bash
snakemake -s Snakefile_LRS --cores 16 --use-conda
```

Again, the workflow should be launched from:

```text
SV-PIPELINE/snakemake_pipelines/lrs
```

because `Snakefile_LRS` loads:

```python
configfile: "config_lrs.yaml"
```

---

# 13. Recommended validation sequence

Do not start a full production analysis immediately after cloning the repository.

Use:

```text
1. Clone repository
        ↓
2. Edit configuration
        ↓
3. Check reference files
        ↓
4. Check BAMs
        ↓
5. Validate Conda environments
        ↓
6. Snakemake --lint
        ↓
7. Snakemake -n
        ↓
8. Test one sample/rule where appropriate
        ↓
9. Run production workflow
```

Useful commands:

```bash
git status
```

```bash
conda --version
```

```bash
snakemake --version
```

```bash
python -m py_compile scripts/*.py
```

```bash
snakemake -s snakemake_pipelines/srs/Snakefile_SRS --lint
```

```bash
snakemake -s snakemake_pipelines/lrs/Snakefile_LRS --lint
```

---

# 14. Important reference resources

The repository contains or expects several different classes of reference data.

## Genome reference

```text
hg38.fa
```

All major workflows are designed around GRCh38/hg38.

## Optic-neuropathy panel

```text
PANEL_OA/optic_neuropathy_panel.bed
PANEL_OA/optic_neuropathy_genes.txt
```

Used for targeted/panel-oriented analysis and candidate-gene filtering.

## HPO reference

```text
reference/optic_neuropathy_hpo_terms.tsv
```

Used for the human gene–phenotype discovery layer.

## STR catalog

```text
reference/variant_catalog_hg38_repeats.json
```

Used by the SRS ExpansionHunter workflow.

## LRS repeat/transposable-element references

LRS additionally requires resources such as:

```text
human_GRCh38_no_alt_analysis_set.trf.bed
teref.ont.human.fa
```

These must exist at the configured locations.

---

# 15. Important tool-specific details

## Sniffles2

The current LRS workflow uses:

```text
--minsupport 3
--no-qc
--allow-overwrite
--min-alignment-length 100
--output-rnames
--tandem-repeats <VNTR BED>
```

Mosaic calling is configurable:

```yaml
sniffles_mosaic: false
```

so mosaic mode is not silently enabled by default.

## cuteSV

The current workflow uses:

```text
--min_support 3
--genotype
```

with insertion/deletion clustering parameters configured in the Snakefile.

## Jasmine

The two LRS SV callsets are merged with:

```text
Sniffles2 + cuteSV → Jasmine
```

and genotype information is requested.

## VEP

VEP is configured for offline cache-based annotation:

```text
--cache
--offline
--dir_cache <VEP_CACHE_DIR>
--assembly GRCh38
```

The full VEP output is retained, while a panel-focused view is generated separately.

## DeepVariant

SRS can use DeepVariant through either:

- a Conda environment, or
- the configured DeepVariant container

depending on:

```yaml
deepvariant_use_conda:
```

The current default is container-based execution.

---

# 16. Output organization

Outputs are organized per sample according to analysis type.

Typical SRS structure:

```text
<output>/<sample>/
├── coverage/
├── sv/
│   ├── manta/
│   ├── delly/
│   ├── merged/
│   └── vep/
├── mei/
│   └── melt/
├── str/
│   └── expansionhunter/
├── snp_deepvariant/
├── phasing/
└── gene_discovery/
```

Typical LRS structure:

```text
<output>/<sample>/
├── coverage/
├── snp_clair3/
├── phasing/
├── phasing_longphase/
├── sv/
│   ├── sniffles2/
│   ├── cutesv/
│   ├── merged/
│   ├── insertions/
│   ├── straglr/
│   └── needlr/
├── sv/annotsv/
├── methylation/
├── mei/
│   └── tldr/
└── gene_discovery/
```

The exact output locations are controlled by `path:` in the corresponding configuration file.

---

# 17. Interpretation strategy

The pipeline should be interpreted in layers rather than treating one output as definitive.

```text
SV discovery
     ↓
cross-caller / within-platform consolidation
     ↓
AnnotSV annotation
     ↓
frequency / annotation evidence
     ↓
panel genes
     +
genome-wide non-panel genes
     ↓
human phenotype associations
     ↓
candidate ranking
     ↓
manual clinical/genomic review
```

A high-ranked gene is therefore a **candidate for investigation**, not automatically a disease-causing gene.

Similarly, an AnnotSV annotation or SV call should not be interpreted as pathogenic solely because it overlaps a disease-associated gene.

---

# 18. SRS vs LRS

| Component | SRS | LRS |
|---|---|---|
| Input | Short-read BAM | Long-read BAM |
| Coverage | mosdepth | mosdepth |
| SV callers | Manta + Delly | Sniffles2 + cuteSV |
| SV merge | SURVIVOR | Jasmine |
| Population/reference SV filtering | — | needLR |
| Small variants | DeepVariant | Clair3 |
| Phasing | WhatsHap | WhatsHap + LongPhase |
| SV annotation | AnnotSV | AnnotSV |
| Transcript annotation | VEP | VEP |
| MEI analysis | MELT | TLDR |
| STR analysis | ExpansionHunter | Straglr |
| Insertion classification | — | RepeatMasker |
| Methylation | — | modkit |
| Gene discovery | Yes | Yes |
| Human HPO layer | Yes | Yes |
| Candidate ranking | Yes | Yes |

The workflows therefore share the same downstream interpretation architecture while retaining sequencing-specific discovery tools.

---

# 19. Current configuration warnings

Before production execution, check every machine-specific path.

In particular, the current repository configuration contains paths that are intended to be changed.

### SRS

Check:

```yaml
path:
ref:
samples:
melt_dir:
melt_transposons:
melt_genes_bed:
vep_cache_dir:
```

The current SRS MELT/VEP entries contain placeholder-style prefixes and should not be assumed to exist on another computer.

### LRS

Check:

```yaml
path:
ref:
samples:
modkitenv:
clair3env:
vntr_bed:
vep_cache_dir:
needlr_backend_files:
tldr_elements:
```

These paths depend on the installation and reference layout of the analysis server.

---

# 20. Reproducibility

For reproducible analyses, record:

```text
Git commit
Snakemake version
Conda version
environment YAML versions
reference genome version
reference database versions
SV caller versions
sample configuration
```

A recommended record is:

```bash
git rev-parse HEAD
snakemake --version
conda --version
```

and a copy of the exact configuration used for the analysis.

---

# 21. Project status

The repository currently contains:

- SRS structural-variant workflow
- LRS structural-variant workflow
- platform-specific SV discovery
- genome-wide AnnotSV
- optic-neuropathy panel-focused annotation
- human gene/HPO discovery
- non-panel candidate extraction
- integrated candidate ranking
- VEP annotation
- MEI/STR/long-read-specific analyses
- Conda environment definitions for pipeline components

The pipeline should be considered **ready for environment-specific validation**, but not automatically assumed to be production-ready on a new computer until all external tools, references, databases, executable paths, and sample paths have been verified.

---

## 22. Citation and software references

When using this workflow for research, cite the individual software packages used by the relevant pipeline components, including:

- Snakemake
- Manta
- Delly
- SURVIVOR
- Sniffles2
- cuteSV
- Jasmine
- needLR
- AnnotSV
- VEP
- Clair3
- DeepVariant
- WhatsHap
- LongPhase
- mosdepth
- RepeatMasker
- Straglr
- ExpansionHunter
- MELT
- modkit
- TLDR
- Monarch Initiative / Monarch Knowledge Graph

The exact versions should be taken from the corresponding `envs/*.yaml` files used for the run.

---

## 23. Repository

**GitHub:**  
https://github.com/kevR27/SV-PIPELINE

**Main workflows:**

```text
snakemake_pipelines/srs/Snakefile_SRS
snakemake_pipelines/lrs/Snakefile_LRS
```

**Configuration:**

```text
snakemake_pipelines/srs/config_srs.yaml
snakemake_pipelines/lrs/config_lrs.yaml
```

**Gene-discovery scripts:**

```text
scripts/extract_annotsv_genes.py
scripts/extract_nonpanel_genes.py
scripts/monarch_human_gene_phenotypes.py
scripts/rank_sv_gene_candidates.py
```

**Panel:**

```text
PANEL_OA/
```

**Reference resources:**

```text
reference/
```
