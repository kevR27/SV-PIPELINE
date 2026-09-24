# SV-PIPELINE

Long- and short-read whole-genome sequencing workflows for genome-wide structural-variant discovery and interpretation, developed for unresolved rare neurological disease with a particular focus on optic neuropathies and mitochondrial/neuromuscular phenotypes.

The repository is designed around a diagnostic question rather than around a single caller:

> Can whole-genome sequencing identify clinically relevant genomic mechanisms that may be missed by an analysis focused mainly on SNVs, small indels, or conventional CNV testing?

The long-read workflow uses Oxford Nanopore whole-genome data and combines several independent SV callers with gene/phenotype annotation and orthogonal evidence. The optic-neuropathy gene panel is used for interpretation, but it is **not used to restrict genome-wide SV discovery**. This preserves the possibility of finding variants in known disease genes as well as new or unexpected candidate genes.

The active long-read workflow is `snakemake_pipelines/lrs/Snakefile_LRS_update`. A separate downstream workflow, `Snakefile_LRS_postprocess`, adds orthogonal evidence and generates interpretation tables and thesis figures without rerunning variant calling.

---

## 1. Biological hypothesis and analysis questions

The working hypothesis is that a proportion of unresolved neurological/optic-neuropathy cases may be explained by genomic alterations that are difficult to capture or interpret with standard small-variant analysis alone. These can include:

- deletions, duplications, inversions, insertions and breakends;
- large or complex rearrangements;
- tandem-repeat expansions;
- mobile-element insertions;
- variants whose interpretation depends on phasing or local genomic context;
- regulatory or methylation changes that may provide additional biological context.

The pipeline therefore asks several linked questions.

**Question 1 — Is there a structural variant in the genome that could explain the phenotype?**  
Three long-read SV callers are used independently and merged into a patient-level master SV callset.

**Question 2 — Is the event technically credible?**  
Caller-specific read support, FILTER status, caller concordance, size, QC flags and sequencing depth are retained rather than reduced to a single binary decision.

**Question 3 — Does the event affect a biologically relevant gene or region?**  
The master callset is annotated genome-wide with AnnotSV and supplemented with VEP.

**Question 4 — Is the event rare enough to be compatible with a rare-disease hypothesis?**  
needLR is used as a supplementary ONT population-frequency source and is matched back to the master SVs.

**Question 5 — Does the affected gene fit the patient's disease biology?**  
Panel membership, OMIM/GenCC information and offline Monarch/HPO phenotype relationships are combined for prioritization.

**Question 6 — Is there additional long-read evidence that changes the interpretation?**  
Straglr, TLDR, LongPhase/WhatsHap and modkit provide repeat, mobile-element, phasing and methylation context.

The result is not an automatic pathogenicity classifier. The pipeline produces an auditable evidence framework that helps select variants for expert review, orthogonal confirmation and clinical interpretation.

---

## 2. Why whole-genome discovery is used

The structural-variant arm is intentionally genome-wide. Restricting SV discovery to an optic-neuropathy BED would make the analysis blind to:

- large rearrangements extending outside a panel interval;
- breakpoints outside coding exons;
- regulatory or intergenic events;
- variants affecting genes not yet included in the panel;
- new candidate genes with a phenotype relationship to the patient.

The optic-neuropathy panel is therefore applied **after discovery** as an interpretation layer. A variant can be labelled as a panel-gene event or a non-panel candidate without changing whether it remains in the master callset.

The LRS design assumes standard ONT whole-genome sequencing. Adaptive sampling is not part of the active workflow.

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
thesis-quality plots and candidate summaries
```

---

## 4. Step-by-step rationale

### 4.1 mosdepth — sequencing-depth QC

**Input:** aligned ONT BAM  
**Output:** per-sample depth distribution and summary files.

mosdepth is used before interpretation because low or uneven coverage can explain missed calls or unstable genotype/support estimates. Coverage is not itself evidence that an SV is pathogenic, but it is essential context when comparing samples or investigating why one caller failed to detect a known event.

Main outputs:

```text
<sample>/coverage/<sample>.mosdepth.global.dist.txt
<sample>/coverage/<sample>.mosdepth.summary.txt
```

**Question answered:** Was the genome sequenced deeply and uniformly enough for the downstream result to be interpretable?

---

### 4.2 Clair3 — complementary SNV/indel arm

Clair3 is not the main discovery engine of this project. It is retained as a complementary small-variant layer, currently targeted to the optic-neuropathy BED.

**Why it is included:** an unresolved patient can still carry relevant SNVs/indels, and small variants are useful for phasing with nearby structural variants.

**Output:** compressed small-variant VCF.

```text
<sample>/snp_clair3/<sample>.vcf.gz
```

**Question answered:** Are there small variants in the disease-focused regions that should be considered together with the structural-variant result?

---

### 4.3 WhatsHap — small-variant phasing and haplotagging

WhatsHap phases the Clair3 calls using the patient's long reads and creates a haplotagged BAM.

**Why it is included:** phase can determine whether multiple variants occur on the same or different haplotypes and can support later locus-specific interpretation.

Main outputs:

```text
<sample>/phasing/<sample>.phased.vcf.gz
<sample>/phasing/<sample>.phased.bam
```

**Question answered:** Which small variants can be assigned to the same haplotype?

---

## 5. Genome-wide structural-variant discovery

No single SV caller is treated as complete. Long-read callers differ in how they cluster read signatures, represent breakpoints and apply internal QC. The workflow therefore uses three complementary callers.

### 5.1 Sniffles2

Sniffles2 is a primary long-read SV caller and reports deletions, duplications, insertions, inversions and breakend-type events from long-read alignments.

The main discovery environment is kept separate from the needLR-compatible Sniffles environment. The repository currently pins the main discovery branch to Sniffles2 2.8.1.

The workflow can optionally expose QC-failed Sniffles candidates with `--qc-output-all` when the controlled `COV_VAR` rescue is enabled. Importantly, this does **not** mean that every QC-failed candidate is accepted.

Current rescue logic:

```text
Sniffles FILTER = COV_VAR
AND SVTYPE in {DEL, DUP}
AND read support >= 2
AND |SVLEN| >= 50 kb
        ↓
retain for downstream review
        ↓
mark EVIDENCE_FLAGS = RESCUED_COV_VAR
```

Other Sniffles QC failures remain excluded.

**Why this exists:** validation against known large deletions showed that caller-version-specific coverage QC can remove a biologically real large SV even when the breakpoint and supporting reads are present. The rescue is therefore narrow and auditable rather than a global `--no-qc` policy.

---

### 5.2 cuteSV

cuteSV provides an independent long-read SV signature/clustering approach.

**Why it is included:** agreement between callers increases technical confidence, while disagreement is also informative because different algorithms can detect different event representations.

The ONT-oriented clustering parameters are defined explicitly in the Snakefile. Parameters for very large deletions should be benchmarked against known-positive samples rather than assumed to be optimal for every dataset.

**Question answered:** Does an independent long-read clustering strategy recover the same event, and if not, where does caller sensitivity differ?

---

### 5.3 DELLY long-read mode

DELLY is run in ONT long-read mode.

**Why it is included:** it provides a third, algorithmically distinct source of breakpoint/SV evidence. A variant that is missed by Sniffles2 or cuteSV may still be retained by DELLY, and vice versa.

**Question answered:** Is the event supported by an additional SV-calling model?

---

## 6. Caller normalization and pre-Jasmine evidence filtering

The raw caller VCFs use different INFO/FORMAT fields. `parse_sv_caller_vcf.py` converts them into a common caller-evidence schema containing coordinates, SV type, size, FILTER status and normalized support.

The normalized evidence is then processed by `filter_sv_evidence.py`.

Current hard filters include:

```text
minimum caller support = 2
minimum |SVLEN| = 50 bp
no global upper SV-length cutoff
```

Additional information such as imprecision, low genotype quality or blacklist overlap can be retained as flags rather than automatically deleting the record.

The filtering script streams rows instead of loading the complete genome-wide table into memory, which is important when Sniffles `--qc-output-all` creates large evidence files.

**Why this step is separate from the caller:** the goal is to make the inclusion criteria visible and comparable across callers before merging.

**Question answered:** Which caller records have enough basic evidence to be presented to the cross-caller merging step?

---

## 7. Jasmine — the master patient SV universe

The filtered caller VCFs are merged per patient with Jasmine using a fixed order:

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

The workflow uses `--allow_intrasample` because all three inputs are different caller representations of the **same biological sample**.

The complete sorted/indexed Jasmine VCF is the **source-of-truth master SV universe**:

```text
<sample>/sv/merged/<sample>_merged_SV.vcf.gz
```

A second VCF containing variants supported by at least the configured number of callers is created as a companion high-confidence evidence set:

```text
<sample>/sv/merged/<sample>_merged_SV.high_confidence.vcf.gz
```

This companion file does **not** replace the complete master callset. A clinically relevant single-caller SV is not automatically false.

**Question answered:** Which caller representations likely describe the same biological event, and how much cross-caller support does each master SV have?

---

## 8. Caller-support summary

`summarize_sv_caller_support.py` decodes Jasmine provenance and generates:

```text
<sample>/sv/merged/<sample>_caller_support_summary.tsv
```

This table records which callers support each master SV and provides the basis for caller-concordance plots.

**Interpretation:** caller count is technical evidence, not a pathogenicity score.

---

## 9. needLR — supplementary population-frequency evidence

needLR is intentionally kept outside the Jasmine master-calling chain.

A dedicated Sniffles2 2.6.2 query VCF is created because the needLR ONT backend was built from a Sniffles2-compatible representation. The cross-caller Jasmine VCF is therefore **not** used as the direct needLR query.

```text
ONT BAM
  ↓
Sniffles2 2.6.2
  ↓
needLR
  ↓
population-frequency evidence
  ↓
coordinate/type matching back to Jasmine master SV
```

Main outputs:

```text
<sample>/sv/needlr/<sample>_needLR_RESULTS.tsv
<sample>/sv/needlr/<sample>_needLR_RESULTS.vcf.gz
```

Important rules:

- needLR is a supplementary annotation source, not a gatekeeper.
- a missing needLR match must **not** be interpreted as allele frequency zero.
- BNDs and SVs >=10 Mb remain in the master analysis even when they are not evaluable by needLR.
- needLR evidence is matched by compatible SV type and genomic coordinates, not merely because two records affect the same gene.

**Question answered:** Is a compatible representation of this SV observed in the available ONT population reference, and at what frequency?

---

## 10. AnnotSV — primary structural-variant annotation

AnnotSV consumes the complete Jasmine master VCF directly.

It provides gene overlap and clinically useful SV annotations and is run genome-wide first. A panel-only view is then derived from the full annotation.

Main outputs:

```text
<sample>/sv/annotsv/<sample>_merged_SV.annotsv.tsv
<sample>/sv/annotsv/<sample>_merged_SV.annotsv.panel_only.tsv
```

**Why this order matters:** panel filtering before annotation would defeat the genome-wide discovery hypothesis.

**Question answered:** Which genes and clinically relevant genomic features are affected by each master SV?

---

## 11. VEP — supplementary consequence annotation

VEP is run on the same complete Jasmine master VCF.

It is supplementary rather than the primary SV annotation source.

Main outputs:

```text
<sample>/sv/vep/<sample>_SV_VEP.txt
<sample>/sv/vep/<sample>_SV_VEP.panel_only.txt
```

**Question answered:** What additional transcript/consequence context can be attached to the master SV representation?

---

## 12. Genome-wide gene and phenotype discovery

After AnnotSV, all genes affected by master SVs are extracted.

The pipeline creates both:

```text
<sample>_SV_genes.txt
<sample>_nonpanel_genes.txt
```

This is intentional. The panel helps recognize known optic-neuropathy genes, but non-panel genes remain available for discovery.

The offline Monarch Knowledge Graph/HPO branch then asks whether affected genes are connected to the disease phenotype anchors.

Outputs include:

```text
<sample>_human_gene_phenotypes.tsv
<sample>_ranked_candidates.tsv
```

The ranking combines biological context such as panel membership and phenotype relationships to prioritize review.

**Important:** `CANDIDATE_CLASS` and phenotype scores are prioritization aids. They are not automatic clinical classifications and should not be treated as substitutes for variant interpretation guidelines.

**Question answered:** Which affected genes are most biologically compatible with the phenotype, including genes outside the original panel?

---

## 13. Integrated SV/gene evidence table

The main interpretation product of the core LRS workflow is:

```text
<sample>/gene_discovery/<sample>_integrated_SV_gene_analysis.tsv
```

The table is defined from the Jasmine master callset. One row is written per **(master SV, overlapping gene)**.

It integrates:

- master SV coordinates, type and length;
- Jasmine caller provenance;
- caller count and `SUPP_VEC`;
- normalized read-support evidence;
- caller evidence/QC flags;
- affected genes;
- AnnotSV matching and annotation;
- needLR allele-frequency evidence and evaluation status;
- OMIM and GenCC information where available;
- optic-neuropathy panel status;
- phenotype score;
- candidate-priority class;
- retained master VCF INFO fields.

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
GENES
NEEDLR_AF / NEEDLR_STATUS
OMIM / GENCC
ANNotsv_Classification
PANEL_STATUS
PHENOTYPE_SCORE
CANDIDATE_CLASS
```

This table is the main bridge between **variant discovery** and **biological interpretation**.

---

## 14. Orthogonal long-read evidence

The following branches are intentionally not counted as additional Jasmine SV callers. They answer different biological questions.

### 14.1 Straglr — tandem-repeat expansions

Straglr searches for tandem-repeat expansions and the results are gene-annotated.

Outputs:

```text
<sample>/sv/straglr/<sample>_straglr.vcf
<sample>/sv/straglr/<sample>_straglr.tsv
<sample>/sv/straglr/<sample>_straglr.annotated.tsv
```

**Question answered:** Is there a repeat-expansion mechanism that a conventional breakpoint-SV callset may not represent adequately?

---

### 14.2 TLDR — mobile-element insertions

TLDR provides an independent mobile-element insertion analysis.

Output:

```text
<sample>/mei/tldr/<sample>.tldr.table.txt
```

**Question answered:** Is a mobile-element insertion present at or near a candidate locus?

---

### 14.3 LongPhase — SNP/SV phasing

LongPhase combines the same-patient small variants, filtered Sniffles SVs, BAM and reference.

Output:

```text
<sample>/phasing_longphase/<sample>.longphase.vcf.gz
```

**Question answered:** Can a candidate SV be placed in haplotypic context with nearby sequence variants?

---

### 14.4 modkit — methylation context

modkit generates CpG methylation evidence from modification-aware ONT BAMs, retaining 5mC and 5hmC as separate modification classes.

Main output:

```text
<sample>/methylation/<sample>.cpg.bedmethyl.gz
```

Methylation is treated as **context**, not as proof that an SV is causal.

**Question answered:** Does a candidate region show an informative local methylation pattern that may support follow-up biological interpretation?

---

## 15. Post-processing: integrating orthogonal and multimodal context

Run `Snakefile_LRS_postprocess` only after the core LRS workflow has produced its integrated table.

The post-processing workflow does not recall variants. It preserves the Jasmine master SV universe and attaches additional evidence.

It creates:

```text
<sample>_integrated_SV_gene_with_orthogonal_evidence.tsv
<sample>_independent_orthogonal_findings.tsv
<sample>_integrated_SV_gene_with_multimodal_context.tsv
<sample>_gene_multimodal_evidence_summary.tsv
```

### Integrated orthogonal-evidence table

Adds coordinate-aware Straglr, TLDR and optional LongPhase evidence to master SVs.

### Independent orthogonal findings

Retains relevant Straglr/TLDR findings that do **not** have a compatible Jasmine counterpart. This prevents biologically interesting repeat or MEI findings from disappearing simply because they are not represented as a conventional master SV.

### Multimodal-context table

Adds nearby WhatsHap-phased small variants and local methylation context. These layers remain distinct from SV caller support.

### Gene-level multimodal summary

Collapses evidence at the gene level for candidate review without double-counting the same master SV simply because it overlaps more than one annotation row.

---

## 16. What information is available at the end?

The pipeline is designed so that the final review can answer, for each candidate:

| Evidence layer | Final information |
| --- | --- |
| Technical QC | sequencing depth and coverage distribution |
| Master SV | chromosome, breakpoints, SV type, SV length |
| Caller provenance | which of Sniffles2/cuteSV/DELLY detected the event |
| Read support | normalized per-caller supporting evidence |
| Caller QC | PASS/FAIL-derived flags and controlled rescue annotations |
| Cross-caller confidence | `SUPP`, `SUPP_VEC`, caller count |
| Gene effect | overlapping/affected genes from AnnotSV |
| Known-disease context | panel membership, OMIM, GenCC |
| Population evidence | needLR AF/status when evaluable |
| Phenotype relevance | Monarch/HPO relationship and phenotype score |
| Candidate priority | discovery-oriented candidate class |
| Repeat evidence | Straglr repeat-expansion findings |
| MEI evidence | TLDR mobile-element findings |
| Haplotype context | WhatsHap/LongPhase phasing |
| Epigenetic context | modkit methylation information |
| Final review level | SV-level, SV-gene-level and gene-level integrated tables |

The intended end point is therefore not simply “a VCF.” It is a traceable chain:

```text
raw long-read evidence
        ↓
caller-specific SV hypotheses
        ↓
transparent QC
        ↓
patient-level master SV
        ↓
gene + clinical annotation
        ↓
population-frequency evidence
        ↓
phenotype prioritization
        ↓
repeat / MEI / phasing / methylation context
        ↓
candidate variants for expert review and confirmation
```

---

## 17. Source-of-truth and interpretation rules

These rules are central to the pipeline design.

1. The complete Jasmine VCF is the master LRS structural-variant universe.
2. A variant is not removed merely because needLR, VEP or another annotation layer cannot evaluate it.
3. Multi-caller support increases technical evidence but does not define pathogenicity.
4. Single-caller variants remain available for interpretation.
5. The high-confidence multi-caller VCF is a companion set, not the master callset.
6. Panel membership is an annotation, not a genome-wide discovery restriction.
7. needLR absence is not equivalent to allele frequency zero.
8. BNDs and very large SVs remain in the master analysis even if a supplementary tool has a restricted evaluation range.
9. Straglr, TLDR, phasing and methylation are orthogonal/contextual layers and do not inflate Jasmine caller count.
10. Candidate ranking is prioritization, not an automatic diagnostic classification.
11. Known-positive samples should be used to benchmark sensitivity before applying parameter changes to the full cohort.

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

Dry run first:

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

The output root is controlled by `path:` in `config_lrs.yaml`.

---

## 19. Running post-processing

Build the downstream interpretation tables:

```bash
snakemake \
  --snakefile Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8
```

Generate all applicable thesis plots explicitly:

```bash
snakemake \
  --snakefile Snakefile_LRS_postprocess \
  --configfile config_lrs.yaml \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 8 \
  all_thesis_plots
```

Plotting is intentionally separated from biological calling so figures can be regenerated without rerunning the expensive WGS analysis.

---

## 20. Thesis plotting layer

The plotting layer generates publication/thesis-oriented PDF, SVG and high-resolution PNG outputs plus source TSV summaries where applicable.

Current plot groups include:

```text
01_qc/
02_caller_concordance/
03_sv_landscape/
04_population_frequency/
05_candidate_prioritization/
06_phenotype/
07_orthogonal/
    ├── straglr/
    └── tldr/
08_phasing/
09_methylation/
10_integrated_evidence/
```

These figures are intended to answer different questions rather than repeat the same information:

- coverage/QC: is the sample technically comparable?
- caller concordance: which calls are shared or caller-specific?
- SV landscape: what types/sizes/chromosomes dominate the callset?
- population frequency: which evaluable SVs are rare/common?
- candidate prioritization: which affected genes deserve review?
- HPO heatmap: which genes have phenotype relationships?
- repeat/MEI plots: are there orthogonal mechanisms outside the breakpoint-SV model?
- phasing: can candidate alleles be placed on haplotypes?
- methylation: is there candidate-region epigenetic context?
- integrated evidence matrix: how many independent evidence layers support each candidate?

See `plots/README.md` for individual plotting commands.

---

## 21. Recommended validation strategy

Before scaling the workflow to an unresolved cohort, validate it against samples with known structural variants.

Recommended checks:

1. confirm the known event is visible in the raw caller outputs;
2. verify whether it passes or fails each caller's internal QC;
3. inspect the normalized caller evidence and filtering decision;
4. verify that Jasmine merges compatible caller representations correctly;
5. inspect `SUPP_VEC`, `IDLIST` and caller provenance;
6. confirm AnnotSV maps the master event to the expected gene(s);
7. verify needLR status without treating an unmatched call as AF=0;
8. confirm the event survives into the integrated SV/gene table;
9. inspect the orthogonal/context layers;
10. only then apply the same settings to unresolved cases.

This validation strategy is particularly important for very large deletions because caller versions and default size/QC parameters can materially change sensitivity.

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

The same conceptual rules are retained:

- the complete merged VCF is the master SRS SV universe;
- multi-caller support is evidence, not an exclusion requirement;
- the panel remains a downstream interpretation layer;
- ExpansionHunter and MELT are independent orthogonal branches;
- needLR is not applied to SRS;
- the shared integrated-table schema reports needLR as not applicable for SRS.

Run with:

```bash
cd /DATA/casadei7/tools/SV-PIPELINE-main/snakemake_pipelines/srs

snakemake \
  --snakefile Snakefile_SRS \
  --use-conda \
  --conda-prefix /home/casadei7/snakemake_envs/envs/ \
  --cores 32 \
  --printshellcmds
```

The SRS branch provides a framework for later LRS-vs-SRS comparison, including assessment of which SV classes, breakpoint representations or difficult genomic regions benefit most from long-read sequencing.

---

## 23. Repository structure

```text
SV-PIPELINE/
├── envs/                       conda environments
├── PANEL_OA/                   optic-neuropathy panel resources
├── reference/                  local reference/annotation resources
├── scripts/                    parsing, filtering and integration scripts
├── plots/                      thesis/post-processing visualization layer
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

A strong candidate does not have to satisfy every layer, but the final tables allow a reviewer to ask:

```text
Is the SV technically supported?
        ↓
Was it detected by one or several callers?
        ↓
Does it affect a known optic-neuropathy gene?
        ├── yes -> established-disease interpretation
        └── no  -> genome-wide candidate discovery remains possible
        ↓
Is the event rare in the available ONT population reference?
        ↓
Does the affected gene fit the phenotype/HPO profile?
        ↓
Is there repeat, MEI, phase or methylation context?
        ↓
Does the total evidence justify manual review,
orthogonal validation and clinical interpretation?
```

That is the central purpose of the pipeline: to move from a large genome-wide SV callset to a transparent, biologically informed shortlist without discarding potentially relevant variants simply because one caller or one annotation resource cannot evaluate them.
