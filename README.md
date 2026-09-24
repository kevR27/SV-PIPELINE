# SV-PIPELINE

This repository contains the workflow I am developing for my Master's project in Pharmaceutical Biotechnology. The main aim is to analyse whole-genome sequencing data from patients with rare neurological diseases, with a particular focus on optic neuropathies, and to investigate whether structural variants or other genomic alterations can help explain cases that remain unresolved after conventional genetic testing.

The project is mainly based on Oxford Nanopore long-read whole-genome sequencing, with a parallel short-read workflow that can later be used for comparison.

The main question behind the pipeline is quite simple:

> **Can whole-genome sequencing reveal clinically relevant variants that may be missed when the analysis is focused mainly on SNVs, small indels or conventional CNVs?**

For this reason, the workflow does not focus only on one type of variant or one list of genes. Structural variants are first searched genome-wide, and the optic-neuropathy gene panel is used later as an interpretation layer.

The active long-read workflow is:

```text
snakemake_pipelines/lrs/Snakefile_LRS_update
```

A second workflow,

```text
snakemake_pipelines/lrs/Snakefile_LRS_postprocess
```

is used afterwards to combine additional evidence and generate the final interpretation tables and thesis plots.

---

## 1. General idea of the project

Patients with optic neuropathies and other rare neurological disorders can remain genetically unresolved even after exome sequencing or targeted analysis.

One possible reason is that the disease-causing alteration is not a simple SNV or small indel. It could instead be:

- a deletion or duplication;
- an inversion;
- a large insertion;
- a complex rearrangement;
- a tandem-repeat expansion;
- a mobile-element insertion;
- a variant whose effect depends on the haplotype;
- or a genomic alteration associated with changes in local methylation.

Long-read sequencing is useful here because individual reads can span much larger genomic regions and can therefore provide information that is difficult to obtain from short reads alone.

The workflow was therefore designed to answer a series of practical questions:

1. **Is there an SV in the genome that could be relevant for the patient's phenotype?**
2. **Is the call technically supported by the sequencing data?**
3. **Do different SV callers identify the same event?**
4. **Does the SV affect a known optic-neuropathy gene or another biologically interesting gene?**
5. **Is the variant rare in available population data?**
6. **Does the affected gene have a phenotype relationship compatible with the patient?**
7. **Is there any additional evidence from repeats, mobile elements, phasing or methylation?**

The pipeline is not intended to automatically decide whether a variant is pathogenic. Its purpose is to organize the different pieces of evidence so that promising candidates can be reviewed more carefully and, when necessary, confirmed with an independent method.

---

## 2. Why the SV analysis is genome-wide

The optic-neuropathy panel is important for interpretation, but I did not want to use it as a hard restriction during SV calling.

A large rearrangement can start outside a gene and still disrupt it. A breakpoint can fall in an intronic or regulatory region. Also, a patient could carry a variant in a gene that is not yet part of the current panel.

For this reason, the main SV callers work genome-wide.

Afterwards, each variant can be classified as affecting:

```text
a known optic-neuropathy/panel gene
or
a non-panel gene that may still be relevant
```

This allows the analysis to keep both a diagnostic component and a discovery component.

The current LRS workflow uses standard ONT whole-genome sequencing rather than adaptive sampling.

---

## 3. Overview of the long-read workflow

```text
ONT WGS BAM
│
├── mosdepth
│     └── coverage QC
│
├── Clair3
│     └── panel SNV/indel calling
│           └── WhatsHap
│                 └── small-variant phasing
│
├── Sniffles2 ─┐
├── cuteSV     ├── filtering/normalization ──> Jasmine ──> MASTER SV VCF
└── DELLY LR   ┘                                      │
                                                      ├── caller support
                                                      ├── AnnotSV
                                                      ├── VEP
                                                      ├── gene annotation
                                                      ├── HPO/Monarch
                                                      └── integrated SV/gene table

ONT BAM ──> Sniffles2 2.6.2 ──> needLR
                                 └── population-frequency information

ONT BAM ──> Straglr ──> tandem-repeat expansions
ONT BAM ──> TLDR ──> mobile-element insertions
ONT BAM ──> modkit ──> methylation
Clair3 + Sniffles2 + BAM ──> LongPhase ──> SNP/SV phasing

                           ↓

                 final integrated evidence
                           ↓
                    thesis plots/tables
```

---

## 4. Quality control and small variants

### mosdepth

The first step is to check sequencing depth with **mosdepth**.

This is important because a missed SV is not always a biological negative result. Sometimes a region simply has insufficient coverage.

The main outputs are:

```text
<sample>/coverage/<sample>.mosdepth.global.dist.txt
<sample>/coverage/<sample>.mosdepth.summary.txt
```

These files are later used to compare coverage between samples and to interpret possible false-negative calls.

### Clair3

**Clair3** is used as a complementary small-variant caller.

The main focus of the project is structural variation, so Clair3 is currently targeted to the optic-neuropathy regions rather than being used as the main genome-wide discovery tool.

Its output can still be useful because a patient may carry relevant SNVs or indels, and these variants can also help with phasing.

```text
<sample>/snp_clair3/<sample>.vcf.gz
```

### WhatsHap

**WhatsHap** phases the small variants using the long reads and produces both a phased VCF and a haplotagged BAM.

```text
<sample>/phasing/<sample>.phased.vcf.gz
<sample>/phasing/<sample>.phased.bam
```

This becomes useful when I want to understand whether different variants are on the same or on different haplotypes.

---

## 5. Structural-variant calling

One of the main choices in this pipeline was not to rely on a single SV caller.

Different callers use different algorithms and can behave differently for the same event, especially for large or complex variants. Using several callers allows me to compare their results instead of assuming that one caller is always correct.

The three long-read callers are:

```text
Sniffles2
cuteSV
DELLY
```

### Sniffles2

Sniffles2 is one of the main callers used for long-read SV discovery.

The main discovery environment is currently kept separate from the Sniffles version used for needLR.

For the discovery branch, the repository currently uses:

```text
Sniffles2 2.8.1
```

During testing with known positive samples, I observed that some large deletions can be detected by Sniffles but removed by its internal coverage-based QC.

For this reason, the workflow contains a controlled rescue for large calls labelled:

```text
COV_VAR
```

The rescue is intentionally conservative:

```text
FILTER = COV_VAR
SVTYPE = DEL or DUP
support >= 2 reads
|SVLEN| >= 50 kb
```

If these conditions are satisfied, the call can continue in the pipeline but is explicitly marked:

```text
RESCUED_COV_VAR
```

This keeps the decision visible instead of silently treating the call as a normal PASS event.

Other Sniffles QC failures are still excluded.

### cuteSV

cuteSV is used as a second independent long-read caller.

Its clustering approach is different from Sniffles2, so it can provide useful confirmation of an event or reveal differences in caller sensitivity.

I am using known-positive samples to test and tune its behaviour, particularly for very large deletions, because default caller parameters are not necessarily optimal for every dataset.

The point is not to force cuteSV to reproduce Sniffles2. Instead, the comparison helps show which events are consistently recovered and which are caller-dependent.

### DELLY

DELLY is used in long-read ONT mode as the third caller.

It gives another independent representation of the SV evidence.

This is particularly useful because a real event can sometimes be missed by one caller but detected by another. In the known-positive samples, this has already been useful for understanding differences between the callers.

---

## 6. Normalizing and filtering the caller results

Sniffles2, cuteSV and DELLY do not report their evidence in exactly the same format.

For this reason, their VCFs are first converted into a common evidence table using:

```text
scripts/parse_sv_caller_vcf.py
```

The table contains information such as:

```text
SV ID
chromosome
start/end
SV type
SV length
FILTER status
read support
genotype-related information
```

The normalized results are then processed by:

```text
scripts/filter_sv_evidence.py
```

The current basic thresholds are:

```text
minimum read support = 2
minimum SV size = 50 bp
no global maximum SV size
```

I chose not to apply an upper size cutoff because large rearrangements are actually one of the variant classes I am interested in.

Some additional features, such as imprecision or low genotype quality, can be kept as flags instead of automatically removing the event.

The filtering script also keeps the decisions visible in the output so that it is possible to understand why a call was retained or rejected.

---

## 7. Jasmine and the master SV callset

After caller-specific filtering, the three VCFs are merged with **Jasmine**.

The input order is fixed:

```text
1. Sniffles2
2. cuteSV
3. DELLY
```

This is important because Jasmine's `SUPP_VEC` is interpreted according to this order.

For example:

```text
100 = Sniffles2 only
010 = cuteSV only
001 = DELLY only
101 = Sniffles2 + DELLY
111 = all three callers
```

The complete Jasmine VCF is considered the **master SV callset for each patient**:

```text
<sample>/sv/merged/<sample>_merged_SV.vcf.gz
```

I do not remove an SV simply because it was detected by only one caller.

A single-caller event can still be real, especially when dealing with difficult genomic regions or large rearrangements. Caller concordance is therefore treated as one piece of technical evidence rather than as a strict definition of truth.

A separate multi-caller VCF is also produced:

```text
<sample>/sv/merged/<sample>_merged_SV.high_confidence.vcf.gz
```

This is useful as a higher-confidence companion set, but it does not replace the complete Jasmine VCF.

---

## 8. Population frequency with needLR

For a rare disease project, one of the most useful questions is whether a candidate SV is common in the population.

For this I use **needLR**.

needLR is kept as a separate branch because its ONT backend was generated using a Sniffles2-compatible representation.

The workflow therefore uses:

```text
ONT BAM
   ↓
Sniffles2 2.6.2
   ↓
needLR
```

The main results are:

```text
<sample>/sv/needlr/<sample>_needLR_RESULTS.tsv
<sample>/sv/needlr/<sample>_needLR_RESULTS.vcf.gz
```

The needLR result is then matched back to the Jasmine master variants by compatible SV type and genomic coordinates.

This is important because the Jasmine SV and the needLR query SV are not necessarily represented with exactly the same ID or breakpoint.

A missing needLR match is **not automatically considered AF = 0**.

It simply means that the variant could not be matched/evaluated in that branch.

Also, variants that are outside needLR's supported range, such as BNDs or very large SVs, are still kept in the main analysis.

---

## 9. Annotating the structural variants

### AnnotSV

**AnnotSV** is the main annotation tool used for the structural variants.

It is run directly on the complete Jasmine VCF.

The genome-wide output is:

```text
<sample>/sv/annotsv/<sample>_merged_SV.annotsv.tsv
```

A panel-only view is then extracted:

```text
<sample>/sv/annotsv/<sample>_merged_SV.annotsv.panel_only.tsv
```

This order is important.

I first want to know what the SV affects genome-wide. Only afterwards do I ask whether one of the affected genes belongs to the optic-neuropathy panel.

### VEP

**VEP** is also run on the Jasmine VCF as a supplementary annotation layer.

```text
<sample>/sv/vep/<sample>_SV_VEP.txt
```

In this pipeline, AnnotSV remains the main structural-variant annotation source, while VEP provides additional transcript/consequence information.

---

## 10. From a list of genes to biologically relevant candidates

Once the SVs are annotated, the pipeline extracts all affected genes.

Two useful lists are produced:

```text
<sample>_SV_genes.txt
<sample>_nonpanel_genes.txt
```

The second file is important because I do not want the analysis to stop at the current panel.

The affected genes are then compared with phenotype information from an offline Monarch/HPO knowledge graph.

This produces:

```text
<sample>_human_gene_phenotypes.tsv
<sample>_ranked_candidates.tsv
```

The aim is to prioritize genes that make biological sense for the phenotype.

For example, a non-panel gene should not automatically be ignored if it has a strong relationship with mitochondrial function, retinal ganglion-cell biology, neurodegeneration or relevant HPO terms.

The ranking is therefore a way to decide **what deserves closer inspection**, not an automatic pathogenicity classification.

---

## 11. Main integrated output

The most important table produced by the core LRS workflow is:

```text
<sample>/gene_discovery/<sample>_integrated_SV_gene_analysis.tsv
```

Each row represents a master SV together with one of the genes it overlaps.

The table combines information from the different parts of the pipeline, including:

```text
SV_ID
chromosome/start/end
SVTYPE
SVLEN
which callers detected the event
caller count
SUPP / SUPP_VEC
read support
caller QC flags
affected gene
AnnotSV information
needLR population frequency/status
OMIM
GenCC
panel status
phenotype score
candidate class
```

This is the table I use to move from:

```text
"there are thousands of SVs in this genome"
```

to:

```text
"these are the events that are worth investigating in more detail"
```

---

## 12. Additional long-read evidence

Not every interesting genomic event is best represented as a classical deletion, duplication, insertion or inversion.

For this reason, I also keep several independent analysis branches.

### Straglr

**Straglr** is used to investigate tandem-repeat expansions.

```text
<sample>/sv/straglr/<sample>_straglr.annotated.tsv
```

This is useful because repeat expansions can be disease-causing but may not appear cleanly in the main breakpoint-based SV callset.

### TLDR

**TLDR** is used for mobile-element insertions.

```text
<sample>/mei/tldr/<sample>.tldr.table.txt
```

This provides another type of variant evidence that is not completely covered by the three main SV callers.

### LongPhase

**LongPhase** uses the small variants, filtered Sniffles SVs and long-read BAM to provide SNP/SV phasing.

```text
<sample>/phasing_longphase/<sample>.longphase.vcf.gz
```

This can help place a structural variant in the same haplotypic context as nearby variants.

### modkit

**modkit** is used to obtain CpG methylation information from modification-aware ONT data.

```text
<sample>/methylation/<sample>.cpg.bedmethyl.gz
```

I consider methylation as additional biological context rather than direct confirmation that an SV is pathogenic.

It can become particularly useful when a candidate region has already been identified and I want to look at the local epigenetic pattern.

---

## 13. Post-processing and multimodal integration

The main calling workflow and the interpretation/plotting workflow are separated on purpose.

After the core LRS analysis is complete, I run:

```text
Snakefile_LRS_postprocess
```

This avoids rerunning the expensive variant-calling steps every time I change a plot or add a new interpretation layer.

The post-processing workflow creates:

```text
<sample>_integrated_SV_gene_with_orthogonal_evidence.tsv
<sample>_independent_orthogonal_findings.tsv
<sample>_integrated_SV_gene_with_multimodal_context.tsv
<sample>_gene_multimodal_evidence_summary.tsv
```

The idea is to keep the Jasmine SV callset as the backbone and then add other information around it.

At the same time, Straglr or TLDR findings that do not match a Jasmine SV are not simply thrown away. They are kept separately as independent findings.

This is useful because a repeat expansion or mobile-element insertion may be biologically important even if it is not represented by the main SV callers in exactly the same way.

---

## 14. What I expect to obtain at the end

For each potentially interesting event, I want to be able to answer questions such as:

| Question | Information available |
| --- | --- |
| Was the sample sufficiently covered? | mosdepth |
| What kind of SV is it? | SVTYPE, SVLEN, breakpoints |
| Which caller found it? | Sniffles2 / cuteSV / DELLY |
| How many callers agree? | CALLER_COUNT, SUPP, SUPP_VEC |
| Is there enough read support? | normalized caller support |
| Did the caller flag any QC issue? | evidence flags |
| Which gene is affected? | AnnotSV |
| Is it a known optic-neuropathy gene? | PANEL_STATUS |
| Is it associated with disease? | OMIM / GenCC |
| Is it rare in ONT population data? | needLR |
| Does the gene fit the phenotype? | Monarch/HPO |
| Is there repeat-expansion evidence? | Straglr |
| Is there an MEI? | TLDR |
| Can the event be phased? | WhatsHap / LongPhase |
| Is there useful methylation context? | modkit |

So the final product is not only a VCF.

The analysis moves through several levels:

```text
sequencing data
      ↓
SV calls
      ↓
caller QC
      ↓
merged patient SVs
      ↓
gene annotation
      ↓
population frequency
      ↓
phenotype relevance
      ↓
orthogonal evidence
      ↓
candidate variants for manual review
```

---

## 15. Important interpretation choices

There are a few rules I try to keep consistent throughout the workflow.

- The complete Jasmine VCF is the main LRS SV callset.
- Multi-caller support is useful, but a single-caller call is not automatically discarded.
- needLR is used as population evidence, not as a filter that defines which SVs exist.
- No needLR match does not mean AF = 0.
- The optic-neuropathy panel is used for interpretation, not for restricting genome-wide SV discovery.
- Straglr, TLDR, phasing and methylation are independent evidence layers and are not counted as additional Jasmine callers.
- Candidate ranking helps prioritize variants but does not replace clinical interpretation.
- Large SVs are intentionally retained rather than removed by a global upper-size filter.
- Known-positive samples are used to test whether changes in caller versions or parameters improve or reduce sensitivity.

---

## 16. Validation with known-positive samples

Before applying the workflow to unresolved patients, I am testing it on samples in which a structural variant is already known.

This is especially useful for large deletions because it allows me to check the complete path:

```text
known biological SV
        ↓
raw caller
        ↓
caller QC
        ↓
our filtering
        ↓
Jasmine
        ↓
AnnotSV
        ↓
integrated table
```

For each positive sample, I check:

1. whether Sniffles2 detects the event;
2. whether cuteSV detects it;
3. whether DELLY detects it;
4. whether internal caller QC removes it;
5. whether our evidence filter retains it;
6. whether Jasmine merges compatible calls correctly;
7. whether the expected gene is annotated;
8. whether the event is still present in the final integrated table.

This validation has already shown why it is important not to treat one caller or one software version as absolute truth.

---

## 17. Running the LRS workflow

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main/snakemake_pipelines/lrs

snakemake \
  --snakefile Snakefile_LRS_update \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --rerun-incomplete \
  --printshellcmds
```

I normally check the DAG first with a dry run:

```bash
snakemake \
  --snakefile Snakefile_LRS_update \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --dry-run \
  --printshellcmds
```

The output directory is controlled by:

```yaml
path: ...
```

inside `config_lrs.yaml`.

---

## 18. Running the post-processing workflow

After the core workflow is complete:

```bash
snakemake \
  --snakefile Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8
```

To generate all available thesis plots:

```bash
snakemake \
  --snakefile Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8 \
  all_thesis_plots
```

The plotting step is separate from variant calling so that figures can be regenerated without rerunning the complete WGS pipeline.

---

## 19. Thesis plots

The plotting scripts are organized by biological question:

```text
<sample>/plots/
├── 01_qc/
├── 02_caller_concordance/
├── 03_sv_landscape/
├── 04_population_frequency/
├── 05_candidate_prioritization/
├── 06_phenotype/
├── 07_orthogonal/
│   ├── straglr/
│   └── tldr/
├── 08_phasing/
├── 09_methylation/
└── 10_integrated_evidence/
```

The aim is to generate figures that help answer different parts of the analysis:

- how good the sequencing data are;
- how much the callers agree;
- which SV types and sizes are present;
- which variants are rare;
- which genes are the most interesting;
- how genes relate to the phenotype;
- whether there are repeat or MEI findings;
- whether candidate variants can be phased;
- whether a candidate region has useful methylation information;
- and how the different evidence layers come together for the final candidates.

More details are available in:

```text
plots/README.md
```

---

## 20. Short-read workflow

The repository also contains a parallel short-read WGS workflow.

The idea is similar, but the tools are adapted to Illumina data:

```text
Illumina BAM
│
├── mosdepth
├── DeepVariant -> WhatsHap
│
├── Manta ─┐
└── DELLY ─┴── filtering -> SURVIVOR -> MASTER SRS SV VCF
                                              │
                                              ├── AnnotSV
                                              ├── VEP
                                              ├── caller support
                                              ├── Monarch/HPO
                                              └── integrated table

BAM -> ExpansionHunter
BAM -> MELTv2
```

The SRS workflow can later be compared with the LRS results to investigate which classes of variants are better recovered with long-read sequencing.

This comparison is particularly interesting for large SVs, difficult breakpoints, repetitive regions and other genomic events that are challenging for short reads.

---

## 21. Repository structure

```text
SV-PIPELINE/
├── envs/                       conda environments
├── PANEL_OA/                   optic-neuropathy panel files
├── reference/                  reference and annotation resources
├── scripts/                    parsing/filtering/integration scripts
├── plots/                      plotting and post-processing scripts
└── snakemake_pipelines/
    ├── lrs/
    │   ├── Snakefile_LRS_update
    │   ├── Snakefile_LRS_postprocess
    │   ├── config_lrs.yaml
    │   └── PIPELINE_LOGIC.md
    └── srs/
```

---

## 22. How I interpret a final candidate

The final question is not simply:

```text
"Was this variant called?"
```

Instead, I try to follow the complete evidence:

```text
Is the SV supported by the reads?
        ↓
Which callers detected it?
        ↓
Does it affect a relevant gene?
        ↓
Is the gene already known for optic neuropathy?
        │
        ├── yes -> evaluate in the known disease context
        │
        └── no  -> keep it as a possible genome-wide candidate
        ↓
Is the variant rare?
        ↓
Does the gene fit the phenotype?
        ↓
Is there useful repeat, MEI, phasing or methylation information?
        ↓
Is the total evidence strong enough to justify
manual review and orthogonal confirmation?
```

This is the main purpose of the pipeline: not simply to generate more variants, but to make the large amount of information from whole-genome long-read sequencing easier to interpret in a rare-disease context.
