# SV-PIPELINE

This repository contains long- and short-read whole-genome sequencing workflows for structural-variant (SV) discovery and interpretation in rare neurological disease, with a particular focus on optic neuropathies and mitochondrial/neuromuscular disorders.

The main question behind the pipeline is:

> **Can whole-genome sequencing identify clinically relevant genomic changes that may be missed by analyses focused mainly on SNVs, small indels, or conventional CNVs?**

The long-read workflow uses Oxford Nanopore whole-genome sequencing and combines several independent SV callers with gene annotation, phenotype information, population-frequency data, and additional long-read evidence.

The optic-neuropathy gene panel is used during interpretation, but it is **not used to restrict genome-wide SV discovery**. This keeps the analysis open to both known disease genes and new or unexpected candidate genes.

The active long-read workflow is:

```text
snakemake_pipelines/lrs/Snakefile_LRS_update
```

A separate downstream workflow,

```text
snakemake_pipelines/lrs/Snakefile_LRS_postprocess
```

adds orthogonal evidence and generates integrated interpretation tables and plots without rerunning the main variant-calling steps.

---

## 1. Biological hypothesis and main questions

The working hypothesis is that some unresolved neurological or optic-neuropathy cases may be caused by genomic changes that are difficult to detect or interpret with standard small-variant analysis alone.

These changes can include:

- deletions;
- duplications;
- inversions;
- insertions;
- breakends and complex rearrangements;
- tandem-repeat expansions;
- mobile-element insertions;
- variants whose interpretation depends on phasing;
- genomic changes associated with local methylation differences.

The pipeline is designed to answer a series of connected questions.

**Question 1 — Is there a structural variant that could explain the phenotype?**  
Three long-read SV callers are used independently and their results are merged into one patient-level master SV callset.

**Question 2 — Is the call technically supported?**  
Read support, FILTER status, caller agreement, SV size, QC flags, and sequencing depth are kept as evidence.

**Question 3 — Does the SV affect a relevant gene or genomic region?**  
The master SV callset is annotated genome-wide with AnnotSV and is also analysed with VEP.

**Question 4 — Is the SV rare enough to be compatible with a rare disease?**  
needLR is used as an additional ONT population-frequency source and is matched back to the master SVs.

**Question 5 — Does the affected gene fit the disease phenotype?**  
Panel membership, OMIM/GenCC information, and Monarch/HPO phenotype relationships are used for prioritization.

**Question 6 — Is there additional long-read evidence that helps interpretation?**  
Straglr, TLDR, LongPhase/WhatsHap, and modkit provide repeat, mobile-element, phasing, and methylation information.

The pipeline does not automatically classify a variant as pathogenic. Its purpose is to organize the available evidence so that the most relevant candidates can be reviewed and confirmed more carefully.

---

## 2. Why the SV analysis is genome-wide

The SV analysis is intentionally genome-wide.

If SV calling were restricted only to the optic-neuropathy panel, the workflow could miss:

- large rearrangements that extend outside a panel region;
- breakpoints outside coding exons;
- regulatory or intergenic variants;
- variants affecting genes that are not yet part of the panel;
- new candidate genes that may still be related to the phenotype.

For this reason, the optic-neuropathy panel is applied **after SV discovery** as an interpretation layer.

A variant can therefore be labelled as affecting a known panel gene or a non-panel gene without being removed from the master callset.

The current LRS workflow uses standard ONT whole-genome sequencing. Adaptive sampling is not part of the active workflow.

---

## 3. Core LRS workflow

```text
ONT WGS BAM
│
├── mosdepth
│     └── sequencing-depth QC
│
├── Clair3
│     └── panel SNV/indel calling
│           └── WhatsHap
│                 └── small-variant phasing + haplotagged BAM
│
├── Sniffles2 ─┐
├── cuteSV     ├── caller normalization/QC ──> Jasmine ──> MASTER SV VCF
└── DELLY LR   ┘                                  │
                                                  ├── caller-support summary
                                                  ├── AnnotSV genome-wide
                                                  │      └── panel-derived view
                                                  ├── VEP supplementary annotation
                                                  ├── gene extraction
                                                  ├── Monarch/HPO prioritization
                                                  └── integrated SV/gene evidence table

ONT BAM ──> dedicated Sniffles2 v2.6.2 ──> needLR
                                             └── population-frequency evidence
                                                 matched back to MASTER SVs

ONT BAM ──> Straglr ──> tandem-repeat evidence
ONT BAM ──> TLDR ──> mobile-element insertion evidence
ONT BAM ──> modkit ──> CpG methylation evidence
Clair3 + filtered Sniffles2 + BAM ──> LongPhase ──> SNP/SV phasing

MASTER integrated table
        +
Straglr / TLDR / LongPhase / WhatsHap / methylation
        ↓
post-processing interpretation tables
        ↓
final plots and candidate summaries
```

---

## 4. Step-by-step rationale

### 4.1 mosdepth — sequencing-depth QC

**Input:** aligned ONT BAM  
**Output:** depth distribution and coverage summary for each sample.

mosdepth is used to check whether sequencing coverage is sufficient and reasonably uniform.

This is important because a missed SV does not always mean that the variant is absent. Low or uneven coverage can also reduce caller sensitivity or make support estimates less reliable.

Main outputs:

```text
<sample>/coverage/<sample>.mosdepth.global.dist.txt
<sample>/coverage/<sample>.mosdepth.summary.txt
```

**Question answered:** Is the sequencing coverage sufficient for the downstream analysis?

---

### 4.2 Clair3 — complementary SNV/indel analysis

Clair3 is used as a complementary small-variant layer and is currently focused on the optic-neuropathy regions.

This is useful because an unresolved case can still contain relevant SNVs or small indels, and these variants can also help with phasing.

Main output:

```text
<sample>/snp_clair3/<sample>.vcf.gz
```

**Question answered:** Are there small variants in the disease-focused regions that should be considered together with the SV results?

---

### 4.3 WhatsHap — small-variant phasing and haplotagging

WhatsHap phases the Clair3 variants using the long reads and also creates a haplotagged BAM.

Main outputs:

```text
<sample>/phasing/<sample>.phased.vcf.gz
<sample>/phasing/<sample>.phased.bam
```

This helps determine whether nearby variants are on the same or on different haplotypes.

**Question answered:** Which small variants can be assigned to the same haplotype?

---

## 5. Genome-wide structural-variant discovery

No single SV caller is expected to detect every true event.

Long-read callers use different algorithms, clustering strategies, and QC rules. Because of this, the same biological SV can be represented differently or missed by one caller.

The workflow therefore uses three complementary callers:

```text
Sniffles2
cuteSV
DELLY
```

### 5.1 Sniffles2

Sniffles2 is one of the main long-read SV callers in the workflow. It detects deletions, duplications, insertions, inversions, and breakend-type events from long-read alignments.

The main discovery environment is kept separate from the Sniffles version used for needLR.

The discovery branch currently uses:

```text
Sniffles2 2.8.1
```

The workflow can also expose Sniffles calls that fail internal QC when the controlled `COV_VAR` rescue is enabled.

This does **not** mean that all failed calls are accepted.

The current rescue logic is:

```text
FILTER = COV_VAR
AND SVTYPE = DEL or DUP
AND read support >= 2
AND |SVLEN| >= 50 kb
        ↓
retain the call
        ↓
mark it as RESCUED_COV_VAR
```

Other Sniffles QC failures remain excluded.

**Why this was added:** testing on known large deletions showed that some true events can be detected by Sniffles but removed by its internal coverage-based QC. The rescue is therefore limited to this specific situation and remains visible in the evidence table.

---

### 5.2 cuteSV

cuteSV is used as a second independent long-read SV caller.

Its clustering strategy is different from Sniffles2, so it provides another view of the same sequencing data.

Agreement between callers can increase technical confidence. Differences between callers are also useful because they show where sensitivity or breakpoint representation changes between algorithms.

The ONT-specific clustering parameters are defined in the Snakefile. Very large deletions are validated using known-positive samples so that the chosen settings can be checked against real expected events.

**Question answered:** Does another long-read caller recover the same event, and how similar is its representation?

---

### 5.3 DELLY long-read mode

DELLY is run in ONT long-read mode and provides a third independent source of SV evidence.

A real event can be missed by one caller but still be detected by another. For this reason, DELLY is used together with Sniffles2 and cuteSV rather than as a replacement for either one.

**Question answered:** Is the same event supported by another SV-calling method?

---

## 6. Caller normalization and pre-Jasmine evidence filtering

The raw VCFs from Sniffles2, cuteSV, and DELLY do not use exactly the same INFO and FORMAT fields.

`parse_sv_caller_vcf.py` converts them into a common evidence table containing information such as:

- coordinates;
- SV type;
- SV length;
- FILTER status;
- normalized read support;
- genotype-related fields.

The normalized evidence is then processed by:

```text
scripts/filter_sv_evidence.py
```

Current basic filters include:

```text
minimum caller support = 2
minimum |SVLEN| = 50 bp
no global upper SV-length cutoff
```

Some features, such as imprecision, low genotype quality, or blacklist overlap, can be kept as flags instead of automatically removing the call.

The filtering step also records why a call was retained or rejected.

The script processes the evidence row by row so that large genome-wide files, including Sniffles `--qc-output-all` results, can be handled without loading the complete table into memory.

**Question answered:** Which caller records have enough basic evidence to enter the merging step?

---

## 7. Jasmine — master patient SV callset

After caller-specific filtering, the three VCFs are merged with Jasmine.

The input order is fixed:

```text
1. Sniffles2
2. cuteSV
3. DELLY
```

This order defines the meaning of Jasmine `SUPP_VEC`.

Examples:

```text
100 = Sniffles2 only
010 = cuteSV only
001 = DELLY only
101 = Sniffles2 + DELLY
111 = all three callers
```

The workflow uses `--allow_intrasample` because all three VCFs come from the same biological sample.

The complete sorted and indexed Jasmine VCF is considered the **master SV callset** for each patient:

```text
<sample>/sv/merged/<sample>_merged_SV.vcf.gz
```

A second VCF containing variants supported by at least the configured number of callers is also generated:

```text
<sample>/sv/merged/<sample>_merged_SV.high_confidence.vcf.gz
```

This second file is a companion high-confidence set. It does not replace the complete Jasmine VCF.

A single-caller SV is therefore not automatically considered false.

**Question answered:** Which caller records most likely represent the same biological event, and how many callers support each master SV?

---

## 8. Caller-support summary

`summarize_sv_caller_support.py` decodes the Jasmine caller information and generates:

```text
<sample>/sv/merged/<sample>_caller_support_summary.tsv
```

This table records which callers support each master SV and is also used for caller-concordance plots.

Caller count is treated as technical evidence. It is not a pathogenicity score.

---

## 9. needLR — population-frequency evidence

needLR is kept separate from the Jasmine calling workflow.

A dedicated Sniffles2 2.6.2 query VCF is generated because the needLR ONT backend was built using a Sniffles2-compatible representation.

The workflow is:

```text
ONT BAM
  ↓
Sniffles2 2.6.2
  ↓
needLR
  ↓
population-frequency evidence
  ↓
match back to Jasmine master SV
```

Main outputs:

```text
<sample>/sv/needlr/<sample>_needLR_RESULTS.tsv
<sample>/sv/needlr/<sample>_needLR_RESULTS.vcf.gz
```

The needLR results are matched back to the Jasmine SVs using compatible SV type and genomic position.

Important rules:

- needLR is a supplementary annotation source, not a filter that defines the master SV callset;
- a missing needLR match does **not** mean allele frequency = 0;
- BNDs and SVs >=10 Mb remain in the master analysis even when needLR cannot evaluate them;
- needLR records are matched by SV type and coordinates, not only by gene name.

**Question answered:** Is a compatible SV present in the ONT population reference, and at what frequency?

---

## 10. AnnotSV — main SV annotation and SV-database evidence

AnnotSV is the main annotation tool used for the structural variants.

It is run directly on the complete Jasmine master VCF.

The main outputs are:

```text
<sample>/sv/annotsv/<sample>_merged_SV.annotsv.tsv
<sample>/sv/annotsv/<sample>_merged_SV.annotsv.panel_only.tsv
```

The first file contains the genome-wide annotation. The panel-only file is derived afterwards.

This order is important because the genome-wide discovery step should remain independent of the known gene panel.

In addition to gene annotation and CNV-oriented ACMG/ClinGen ranking, the integrated table now exposes AnnotSV's SV-specific benign/pathogenic overlap fields for:

```text
DEL  -> loss evidence
DUP  -> gain evidence
INS  -> insertion evidence
INV  -> inversion evidence
```

The source fields can contain evidence derived from SV resources such as ClinVar, dbVar, ClinGen, gnomAD, DGV and 1000 Genomes, depending on the installed AnnotSV annotation release.

These are kept as separate columns because the ACMG/ClinGen quantitative framework mainly applies to copy-number loss/gain. For insertions and inversions, the database-overlap evidence is therefore especially useful.

The workflow records both the original source text and simple database-overlap flags. An overlap is treated as evidence, not as proof that the patient SV is exactly the same event as the database SV.

**Question answered:** Which genes and genomic features are affected, and is there known benign or pathogenic SV evidence overlapping the event?

---

## 11. VEP — supplementary transcript/consequence annotation

VEP is also run on the complete Jasmine master VCF.

It is used as a supplementary transcript/consequence layer rather than as the main tool for deciding which genes are affected by an SV.

The workflow uses `--flag_pick`, not `--pick`.

`--flag_pick` keeps all transcript/gene consequences and marks VEP's preferred consequence with the `PICK` flag. This is important for large SVs because one SV can span several genes and transcripts. Using `--pick` would keep only one selected consequence and would under-represent the full overlap.

Main outputs:

```text
<sample>/sv/vep/<sample>_SV_VEP.txt
<sample>/sv/vep/<sample>_SV_VEP.panel_only.txt
```

AnnotSV remains the main genome-wide gene-mapping layer. VEP is used to add transcript-level consequence information.

**Question answered:** What transcript-level consequences are associated with the SV, while keeping all affected genes/transcripts visible?

---

## 12. Genome-wide gene and phenotype discovery

After AnnotSV, all genes affected by the master SVs are extracted.

The workflow creates:

```text
<sample>_SV_genes.txt
<sample>_nonpanel_genes.txt
```

The non-panel list is kept because the analysis should not stop at genes that are already known.

The affected genes are then compared with phenotype information from an offline Monarch/HPO knowledge graph.

This produces:

```text
<sample>_human_gene_phenotypes.tsv
<sample>_ranked_candidates.tsv
```

### Gene ranking

The ranking now keeps three types of information separate:

```text
phenotype relevance
        +
curated gene-disease evidence
        +
a small SV-evidence component
```

The phenotype component has the largest weight because the main purpose is to identify genes that fit the neurological/optic-neuropathy phenotype.

The gene-disease component uses GenCC classifications when available. These are kept as readable terms such as:

```text
Definitive
Strong
Moderate
Supportive
Limited
Animal Model Only
Disputed
Refuted
No known disease relationship
```

A small internal numerical value is attached to these terms only to help order candidates. It is called:

```text
gene_disease_evidence_score
```

This is **not a probability of pathogenicity**.

If GenCC evidence is not available but an OMIM disease relationship is present, a smaller supportive value is used for research prioritization.

The current integrated discovery score is:

```text
phenotype component       0-13
gene-disease evidence     0-4
SV evidence               0-2
                          ----
maximum                    19
```

SV count has only a small contribution because several SVs affecting the same gene do not automatically make that gene more likely to be disease-causing.

Panel membership remains a separate category and is not added directly to the numerical score. This avoids automatically forcing known panel genes above potentially relevant non-panel genes.

`CANDIDATE_CLASS`, `gene_disease_evidence_score`, and the integrated discovery score are prioritization tools. They are not automatic clinical classifications.

**Question answered:** Which affected genes are most compatible with the phenotype and with known human disease evidence, including genes outside the original panel?

---

## 13. Integrated SV/gene evidence table

The main interpretation output of the core LRS workflow is:

```text
<sample>/gene_discovery/<sample>_integrated_SV_gene_analysis.tsv
```

Each row represents one master SV together with one overlapping gene.

The table combines:

- master SV coordinates, type, and length;
- caller provenance;
- caller count and `SUPP_VEC`;
- normalized read support;
- caller QC/evidence flags;
- affected genes;
- needLR population-frequency evidence;
- OMIM and GenCC disease information;
- gene-disease evidence level and score;
- ClinGen haploinsufficiency (HI) and triplosensitivity (TS);
- AnnotSV ranking score and ranking criteria;
- ACMG CNV class for deletions and duplications;
- panel status;
- phenotype and discovery scores;
- relevant INFO fields from the master VCF.

Representative columns include:

```text
SV_ID
CHROM / START / END
CHR2 / POS2
SVTYPE / SVLEN

CALLERS / CALLER_COUNT
SUPP / SUPP_VEC
CALLER_READ_SUPPORT
CALLER_EVIDENCE_FLAGS

NEEDLR_AF / NEEDLR_STATUS

GENES
OMIM
GENCC
GENCC_DISEASE
GENCC_MOI
GENE_DISEASE_EVIDENCE_SCORE
GENE_DISEASE_EVIDENCE_LEVEL
GENE_DISEASE_EVIDENCE_CONFLICT

CLINGEN_HI
CLINGEN_TS
DOSAGE_RELEVANCE

SV_DATABASE_EVIDENCE_SCOPE
SV_PATHOGENIC_DB_SOURCE
SV_PATHOGENIC_DB_COORD
SV_PATHOGENIC_DB_PHENOTYPE
SV_PATHOGENIC_DB_HPO
SV_BENIGN_DB_SOURCE
SV_BENIGN_DB_COORD
SV_BENIGN_DB_AFMAX
SV_DB_CLINVAR_OVERLAP
SV_DB_DBVAR_OVERLAP
SV_DB_GNOMAD_OVERLAP
SV_DB_DGV_OVERLAP
SV_DB_1000G_OVERLAP
SV_DB_CLINGEN_OVERLAP

ANNOTSV_RANKING_SCORE
ANNOTSV_RANKING_CRITERIA
ACMG_CNV_SCORE
ACMG_CNV_CLASS

PANEL_STATUS
PHENOTYPE_SCORE
SV_EVIDENCE_SCORE
INTEGRATED_DISCOVERY_SCORE
CANDIDATE_CLASS
```

For deletions, the dosage interpretation uses the ClinGen **HI** score. For duplications, it uses the ClinGen **TS** score.

The ACMG CNV fields are only treated as formal CNV pathogenicity evidence for `DEL` and `DUP` events. Other SV types keep the AnnotSV ranking information but are not forced into a CNV-specific interpretation framework.

This table connects the SV discovery results with gene relevance, dosage sensitivity, and variant-level pathogenicity evidence.

---

## 14. Orthogonal long-read evidence

The following analyses are kept separate from the main Jasmine caller count because they answer different biological questions.

### 14.1 Straglr — tandem-repeat expansions

Straglr is used to search for tandem-repeat expansions and the results are then linked to genes.

Main outputs:

```text
<sample>/sv/straglr/<sample>_straglr.vcf
<sample>/sv/straglr/<sample>_straglr.tsv
<sample>/sv/straglr/<sample>_straglr.annotated.tsv
```

**Question answered:** Is there a repeat-expansion event that may not be represented well by the main breakpoint-based SV callers?

---

### 14.2 TLDR — mobile-element insertions

TLDR is used to identify mobile-element insertions.

Main output:

```text
<sample>/mei/tldr/<sample>.tldr.table.txt
```

**Question answered:** Is there a mobile-element insertion at or near a candidate locus?

---

### 14.3 LongPhase — SNP/SV phasing

LongPhase combines the small variants, filtered Sniffles SVs, the BAM file, and the reference genome.

Main output:

```text
<sample>/phasing_longphase/<sample>.longphase.vcf.gz
```

**Question answered:** Can a candidate SV be placed on the same haplotype as nearby sequence variants?

---

### 14.4 modkit — methylation context

modkit is used to obtain CpG methylation information from modification-aware ONT BAMs.

Main output:

```text
<sample>/methylation/<sample>.cpg.bedmethyl.gz
```

Methylation is treated as additional biological context, not as direct proof that an SV is disease-causing.

**Question answered:** Does the candidate region show a methylation pattern that may be useful for follow-up interpretation?

---

## 15. Post-processing and multimodal integration

The post-processing workflow is run after the main LRS workflow is complete.

It does not call new master SVs. Instead, it keeps the Jasmine master callset as the backbone and adds other evidence around it.

It creates:

```text
<sample>_integrated_SV_gene_with_orthogonal_evidence.tsv
<sample>_independent_orthogonal_findings.tsv
<sample>_integrated_SV_gene_with_multimodal_context.tsv
<sample>_gene_multimodal_evidence_summary.tsv
```

### Integrated orthogonal-evidence table

Adds coordinate-aware Straglr, TLDR, and optional LongPhase evidence to the master SVs.

### Independent orthogonal findings

Keeps relevant Straglr or TLDR findings that do not have a compatible Jasmine SV.

This is useful because a repeat expansion or mobile-element insertion may still be relevant even if it is not represented by the main SV callers in exactly the same way.

### Multimodal-context table

Adds nearby WhatsHap-phased small variants and local methylation information.

These data are kept separate from SV caller support.

### Gene-level multimodal summary

Combines the different evidence layers at gene level while avoiding repeated counting of the same master SV.

---

## 16. What information is available at the end?

For each candidate, the final analysis can provide:

| Evidence layer | Final information |
| --- | --- |
| Technical QC | sequencing depth and coverage distribution |
| Master SV | chromosome, breakpoints, SV type, SV length |
| Caller provenance | which of Sniffles2/cuteSV/DELLY detected the event |
| Read support | normalized per-caller support |
| Caller QC | QC flags and controlled rescue information |
| Caller agreement | `SUPP`, `SUPP_VEC`, caller count |
| Gene effect | genes affected by the SV |
| Known-disease context | panel membership, OMIM, GenCC |
| Population evidence | needLR AF/status when available |
| Phenotype relevance | Monarch/HPO relationship and phenotype score |
| Candidate priority | discovery-oriented candidate class |
| Repeat evidence | Straglr findings |
| MEI evidence | TLDR findings |
| Haplotype context | WhatsHap/LongPhase |
| Epigenetic context | modkit methylation information |
| Final interpretation | SV-level, SV-gene-level, and gene-level integrated tables |

The final result is therefore more than a VCF.

The analysis follows this general path:

```text
raw long-read evidence
        ↓
caller-specific SV calls
        ↓
transparent QC
        ↓
patient-level master SVs
        ↓
gene and clinical annotation
        ↓
population-frequency evidence
        ↓
phenotype prioritization
        ↓
repeat / MEI / phasing / methylation context
        ↓
candidate variants for review and confirmation
```

---

## 17. Main interpretation rules

The following rules are central to the workflow:

1. The complete Jasmine VCF is the master LRS SV callset.
2. A variant is not removed simply because needLR, VEP, or another annotation tool cannot evaluate it.
3. Support from several callers increases technical confidence but does not define pathogenicity.
4. Single-caller variants remain available for interpretation.
5. The high-confidence multi-caller VCF is a companion file, not the master callset.
6. Panel membership is used for interpretation, not as a restriction during genome-wide SV discovery.
7. No needLR match does not mean allele frequency = 0.
8. BNDs and very large SVs remain in the master analysis even when a supplementary tool has size or representation limits.
9. Straglr, TLDR, phasing, and methylation are additional evidence layers and do not increase the Jasmine caller count.
10. Candidate ranking is used for prioritization, not as an automatic diagnostic classification.
11. The gene-disease evidence score is an internal ranking value derived from curated GenCC/OMIM evidence; it is not a probability that a gene is pathogenic.
12. For deletions and duplications, AnnotSV/ACMG CNV score and class are kept separately from the gene-priority score.
13. ClinGen HI and TS are interpreted according to SV direction: HI for deletions and TS for duplications.
14. SV-specific database overlaps from AnnotSV are retained separately from ACMG CNV classification, especially for insertions and inversions.
15. VEP retains all transcript/gene consequences with `--flag_pick`; AnnotSV remains the main source for genome-wide SV gene mapping.
16. Known-positive samples are used to test sensitivity before parameter changes are applied to the wider dataset.

---

## 18. Running the active LRS workflow

From the LRS workflow directory:

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

A dry run can be used first to inspect the workflow:

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

The output directory is controlled by the `path:` setting in `config_lrs.yaml`.

---

## 19. Running the post-processing workflow

Build the downstream interpretation tables with:

```bash
snakemake \
  --snakefile Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8
```

Generate all available plots with:

```bash
snakemake \
  --snakefile Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8 \
  all_thesis_plots
```

The plotting step is kept separate from the biological analysis so that figures can be regenerated without rerunning the complete WGS workflow.

---

## 20. Plotting layer

The plotting scripts are organized by analysis type:

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

These plots are used to summarize:

- sequencing coverage;
- caller agreement;
- SV type and size distribution;
- population frequency;
- candidate genes;
- phenotype relationships;
- repeat and MEI findings;
- phasing;
- methylation;
- integrated evidence across candidates.

More details are available in:

```text
plots/README.md
```

---

## 21. Recommended validation strategy

Before applying the workflow to unresolved cases, it should be tested using samples with known structural variants.

For each known event, the following should be checked:

1. whether the event is present in the raw caller output;
2. whether it passes or fails the caller's internal QC;
3. whether the normalized evidence is interpreted correctly;
4. whether the pre-Jasmine filter retains the event;
5. whether Jasmine merges compatible caller representations correctly;
6. whether `SUPP_VEC`, `IDLIST`, and caller provenance are correct;
7. whether AnnotSV identifies the expected gene or region;
8. whether needLR status is interpreted correctly;
9. whether the event remains in the integrated table;
10. whether the additional evidence layers provide useful context.

This step is especially important for very large deletions because caller versions and default parameters can change sensitivity.

---

## 22. Short-read workflow

The repository also contains a parallel short-read WGS workflow:

```text
Illumina BAM
│
├── mosdepth
├── DeepVariant -> WhatsHap
│
├── Manta ─┐
└── DELLY ─┴── caller evidence filter -> SURVIVOR -> MASTER SRS SV VCF
                                                   │
                                                   ├── caller-support summary
                                                   ├── AnnotSV
                                                   ├── VEP
                                                   ├── Monarch/HPO
                                                   └── integrated SV/gene table

BAM -> ExpansionHunter
BAM -> MELTv2
```

The same general rules are used:

- the complete merged VCF is the master SRS SV callset;
- multi-caller support is evidence, not an exclusion requirement;
- the gene panel is applied during interpretation;
- ExpansionHunter and MELT remain independent evidence branches;
- needLR is not used for SRS;
- the integrated-table structure remains as similar as possible between LRS and SRS.

The SRS workflow can later be compared with the LRS results to evaluate which types of SVs and genomic regions benefit most from long-read sequencing.

---

## 23. Repository structure

```text
SV-PIPELINE/
├── envs/                       conda environments
├── PANEL_OA/                   optic-neuropathy panel resources
├── reference/                  reference and annotation resources
├── scripts/                    parsing, filtering, and integration scripts
├── plots/                      downstream plotting and interpretation scripts
└── snakemake_pipelines/
    ├── lrs/
    │   ├── Snakefile_LRS_update
    │   ├── Snakefile_LRS_postprocess
    │   ├── config_lrs.yaml
    │   └── PIPELINE_LOGIC.md
    └── srs/
```

---

## 24. Practical interpretation of a final candidate

A useful candidate does not need to be supported by every analysis layer.

The final tables make it possible to ask:

```text
Is the SV technically supported?
        ↓
Was it detected by one or several callers?
        ↓
Does it affect a known optic-neuropathy gene?
        ├── yes -> interpret in the known disease context
        └── no  -> keep it as a possible genome-wide candidate
        ↓
Is the event rare in the available ONT population reference?
        ↓
Does the affected gene fit the phenotype/HPO profile?
        ↓
Is there repeat, MEI, phasing, or methylation information?
        ↓
Does the combined evidence support further review
and orthogonal confirmation?
```

The main purpose of the pipeline is to move from a large genome-wide SV callset to a smaller and more informative set of candidates, while keeping the analysis transparent and without discarding potentially relevant variants simply because one caller or one annotation resource cannot evaluate them.
