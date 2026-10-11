# Candidate-ranking interpretation

This pipeline preserves **research prioritization** and adds a separate clinical
candidate review order. Neither ranking is a diagnostic probability.

The ranking is deliberately divided into evidence axes so that a small numerical
difference in phenotype relevance cannot automatically outrank a much stronger
or weaker structural mechanism.

## Current ranking model

The current LRS ranking model is:

```text
gene relevance tier
    ↓
SV-gene + inheritance/mechanism compatibility
    ↓
technical evidence
    ↓
population-frequency evidence
    ↓
constraint / cohort recurrence / continuous relevance tie-breaks
```

Panel membership is **not** a score. Panel and non-panel genes/events receive
separate within-group ranks.

The main final candidate outputs are:

```text
<sample>_gene_candidates.ranked.tsv
<sample>_sv_gene_candidates.ranked.tsv
```

The authoritative integrated evidence table remains:

```text
<sample>_integrated_SV_gene_analysis.final.tsv.gz
```

## Joint clinical candidate review

The new clinical table asks two questions together: **what classification was
reported for this CNV, and does this actual SV–gene finding merit HON review?**
The answers remain in separate columns. The classification of a large CNV is
an event-level assessment; it does not classify every overlapping gene as the
cause of HON. Each gene's transcript effect, disease mechanism and phenotype
relevance must be reviewed separately.

Outputs inside each sample folder:

| Platform | Clinical table | Existing ranking retained |
| --- | --- | --- |
| LRS | `gene_discovery/final/<sample>_clinical_candidates.ranked.tsv` | Both final gene and SV–gene research tables |
| SRS | `diagnostic_review/<sample>.clinical_candidates.tsv` | The biological-focus `diagnostic_candidates.tsv` table |

Both platforms reuse `build_diagnostic_candidates.py`; LRS calls its shared
function from `finalize_candidate_ranking.py`. The clinical input is the
complete integrated evidence, with final LRS HPO/depth/phase context joined by
exact `SV_ID` + `GENE`. No row is removed, including genes outside the panel,
VUS, unclassified findings, unannotated events and pathogenic findings with
weak or unrelated HON relevance. SRS's optional top-N limit affects its old
review view, not this complete clinical table. Counts of table rows are
SV–gene associations, not counts of unique SVs.

### How the joint order works

P/LP means pathogenic or likely pathogenic **as reported by AnnotSV**, subject
to the audit below. Query matching is a research phenotype comparison, not a
clinical observation or a disease diagnosis.

| Order | `CLINICAL_REVIEW_GROUP` | Meaning |
| --- | --- | --- |
| 1 | `P_LP_HON_MECHANISM_SUPPORTED` | Consistent event classification; strong HON query match; compatible effect and QC genotype; read support; a source-backed disease model with Moderate/Strong/Definitive validity and a compatible mechanism |
| 2 | `P_LP_HON_INCOMPLETE_EVIDENCE` | P/LP with a supporting/strong HON query match, but one or more of those requirements remain unresolved |
| 3 | `HON_RELEVANT_VUS_OR_UNCLASSIFIED` | Supporting/strong HON query match and uncertain/unclassified finding; retained for review without upgrading its class |
| 4 | `P_LP_OTHER_OR_UNRESOLVED_HON` | P/LP with limited or unavailable HON matching; potentially important for another diagnosis |
| 5 | `OTHER_REVIEW` | Other findings, still retained |
| 6 | `CONFLICT_OR_TECHNICAL_REVIEW` | Classification, genotype, mechanism, phenotype, frequency or technical concerns that require resolution; inspect this group separately even if its rows rank lower |
| 7 | `BENIGN_OR_LIKELY_BENIGN` | Consistent reported benign/likely benign class; a strong HPO match does not override it |

Conflict/technical flags route a finding to group 6 before ordinary
prioritisation. A complete, consistent benign class stays in group 7 with any
other concerns still visible. Within a group, the order uses broad query-match
tier, reported class, mechanism support, genotype compatibility and curated
disease-evidence tier, then SV/gene identifiers for deterministic ties. It does
not add points. Panel membership, MitoCarta membership, the old research score,
the magnitude of AnnotSV points within a class, and the number of criteria do
not give extra clinical-ranking bonuses.

For this project, normalized query similarity ≥0.50 is a strong match,
0.25–<0.50 is supporting, and <0.25 is limited. These are transparent research
thresholds, not validated diagnostic cutoffs. Available patient-HPO query
results take precedence, including zero matches; otherwise the shared core
HON query is used, with the existing combined HON score as a legacy fallback.
`CLINICAL_HON_QUERY_SCOPE` records the choice. The cohort query does not imply
that every patient has every queried finding. Exomiser remains a separate
phenotype comparison; agreement is not a required intersection filter.

Gene-level GenCC/OMIM associations and ClinGen dosage sensitivity are useful,
but do not by themselves identify the HON disease model. Disease identifiers,
source, MOI and mechanism from the optional curated allele disease models are
retained when supplied. Without them, `DISEASE_SPECIFIC_HON_MODEL_NOT_AVAILABLE`
is explicit and group 1 is unavailable. Even with a model, the gene-based HPO
match still requires disease-specific clinical review. A heterozygous AR hit
requires another allele; a possible trans partner remains a review candidate.
Only the QC-qualified `ALLELE_GENOTYPE` is used; missing genotype provenance
prevents promotion to group 1. Missing parents do not establish
de novo status or segregation. Balanced/other SVs can merit HON review, but do
not receive a formal DEL/DUP class in this table.

### AnnotSV and the ACMG/ClinGen loss/gain audit

The [2020 ACMG/ClinGen constitutional CNV standard](https://doi.org/10.1038/s41436-019-0686-8)
provides separate copy-number loss and gain frameworks. AnnotSV implements an
[automated adaptation](https://github.com/lgmgeo/AnnotSV/blob/master/share/doc/AnnotSV/ranking.md),
including adaptations for recessive findings; consult the scoring workbook
matching the installed AnnotSV version. This pipeline audits its result; it
does not implement another CNV score to add on top.

`ACMG_CNV_CLASS`, `ACMG_CNV_SCORE`, `ANNOTSV_ACMG_CLASS_RAW`,
`ANNOTSV_RANKING_SCORE` and `ANNOTSV_RANKING_CRITERIA` remain visible. The class
comes from the score and criteria, so these represent one body of evidence.
The audit checks class/point consistency: ≥0.99 P, ≥0.90 LP, >−0.90 to <0.90
VUS, ≤−0.90 LB, and ≤−0.99 B. **0.99 points does not mean a 99% diagnostic
probability.** The reported class is never replaced by a class calculated here.
The audit also checks that the class came from an AnnotSV full/event row and
flags disagreeing event records. A gene split-row class is not silently treated
as the whole-CNV classification. Older tables without that provenance must be
rebuilt or reviewed; no variant callers need rerunning just to rebuild it.

`CLINICAL_CNV_FRAMEWORK_AUDIT` groups the reported criterion codes by section.
No reported code means unavailable/unreported evidence, not proof that the
criterion was evaluated and excluded. `CLINICAL_CNV_MANUAL_OR_UNAVAILABLE`
identifies the following review work:

| CNV section | Existing evidence to inspect | What still needs review |
| --- | --- | --- |
| 1: genomic content | Breakpoints, AnnotSV transcripts/exons/CDS and structural-VEP context | True functional effect, transcript choice and regulatory effects |
| 2: dosage/region evidence | Loss: ClinGen HI; gain: ClinGen TS; pathogenic/benign region matches | Exact allele/region identity, whole versus partial overlap and gain architecture; AD inheritance alone does not establish triplosensitivity |
| 3: gene content | AnnotSV event-level criteria and affected genes | Event-wide coding gene count; do not count repeated gene rows as independent evidence |
| 4: reported cases and population | AnnotSV criteria, available ClinVar/other SV region records and AF sources | Same-allele case/control evidence, phenotype specificity and independent publications; overlapping records can describe different alleles |
| 5: patient/family evidence | Available reliable genotype and source-backed family assessment | Callable family genotypes, segregation, penetrance and clinically observed findings; missing family data stay unknown |

`CLINICAL_PHENOTYPE_IN_CNV_CRITERIA` marks reported case/phenotype criteria that
may already contribute to the class. Their patients' phenotypes are not
automatically this patient's findings. HON matching is used as a separate
relevance axis, never as another ACMG point bonus or a reason to upgrade a VUS.
ClinGen, GenCC/OMIM and ClinVar fields are retained for context/audit without
summing repeated database descriptions of the same underlying evidence.
Mixed pathogenic/benign region records are flagged as context needing allele
identity review, not automatically declared contradictory classifications of
the patient's exact SV. Automated criteria still need an expert audit.

### Flags and practical limits

`CLINICAL_REVIEW_FLAGS` identifies evidence that must be resolved;
`CLINICAL_MISSING_OR_CONTEXT_NOTES` identifies unknowns and interpretation
limits. A technically supported call means computational read support, not
independent validation. Low/ambiguous support, failed caller QC, incompatible
genotype, unclear expected copy change or disease-mechanism conflict cannot be
overridden by a strong HPO match.

The clinical frequency view keeps measured AF ≤1% and unknown AF distinct.
Missing AF is not proof of novelty. needLR coordinate matches remain
provisional; database region AF is not transferred as an exact-allele AF.
These checks do not alter the existing combined rare-or-missing filtering
group, the research ranking or XCI behavior.

This update provides AnnotSV-backed candidate prioritisation. It is not full
Exomiser variant analysis, a new ACMG classification, a calibrated diagnostic
probability or clinical validation. No patient analysis is implied by the code
checks. The OPA1-positive samples 21930 and 25428 are limited controls and do
not establish overall diagnostic performance.

## 1. Phenotype and gene relevance

The old anchor-count score (`min(10, anchor_count * 5)`) is no longer the
primary gene-relevance model.

Generic hereditary optic-neuropathy relevance is calculated with
**directional Resnik phenotype-query coverage** against the broader human
gene-HPO background. Each HON query term is matched to its best gene-associated
HPO term and the result is normalized to the HON query's own self-similarity.
This avoids penalizing pleiotropic genes simply because they have many additional
phenotype annotations. The reference is now split into two roles:

- **core HON ocular terms**, which remain the primary phenotype signal;
- **mitochondrial/syndromic context terms** such as hearing impairment, ataxia,
  lactic acidosis, external ophthalmoplegia, ptosis and exercise intolerance.

The syndromic component is deliberately capped: generic context-only similarity
cannot by itself reach the MODERATE/HIGH phenotype tiers. This prevents broad
mitochondrial features from outranking true optic-neuropathy similarity while
still allowing them to support discovery.

Important fields include:

```text
HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED
HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED
HON_SEMANTIC_SIMILARITY_NORMALIZED
HON_SEMANTIC_METHOD
GENE_RELEVANCE_TIER
GENE_RELEVANCE_DISPLAY_SCORE
GENE_RELEVANCE_PHENOTYPE_SCOPE
```

The display score is balanced between phenotype and curated disease evidence
and is retained for interpretation/tie-breaking. The current generic semantic
method is recorded as:

```text
ASYMMETRIC_RESNIK_QUERY_COVERAGE_NORMALIZED_TO_QUERY_SELF
```

The primary ordering uses the broad relevance tier:

```text
HIGH
MODERATE
SUPPORTING
LIMITED
```

This prevents a 0.5-point difference from dominating mechanism or technical
evidence.

### Patient-specific HPO

If a patient HPO file is configured, the final thesis-analysis layer uses the
same directional Resnik query-coverage logic, with the patient HPO terms as the
query:

```text
HPO_QUERY_RESNIK_NORMALIZED
HPO_BMA_RESNIK_NORMALIZED
FINAL_GENE_RELEVANCE_TIER
FINAL_GENE_RELEVANCE_DISPLAY_SCORE
FINAL_PHENOTYPE_RANKING_SCOPE
```

If patient-specific HPO is unavailable, the pipeline explicitly reports:

```text
GENERIC_HON_FALLBACK_NO_PATIENT_HPO
```

It does not substitute unrecorded syndromic features such as hearing loss,
ataxia, Leigh syndrome or lactic acidosis as though they were present in the
patient.

## 2. Curated gene-disease evidence

GenCC and OMIM evidence remain separate from phenotype similarity.

GenCC strength is first recorded as curated evidence, then reduced when the gene
has little HON phenotype similarity. A positive GenCC record plus contradictory
evidence (for example a Refuted/Disputed record) receives a conflict penalty.

Important fields include:

```text
GENE_DISEASE_EVIDENCE_RAW_SCORE
GENE_DISEASE_EVIDENCE_SCORE
GENE_DISEASE_EVIDENCE_CONFLICT
GENE_DISEASE_HON_CONTEXT_FACTOR
GENE_DISEASE_CONTEXT_SCOPE
```

The current disease-context adjustment is a **gene-level HON semantic proxy**.
It is not yet a disease-ID-specific GenCC-to-HPO join for every record, because
the available integrated rows do not consistently provide a stable disease
identifier that can always be joined to Monarch. The scope field makes this
limitation explicit.

## 3. Dominant/recessive and SV-mechanism compatibility

Inheritance is evaluated at the **SV-gene event** level.

The pipeline uses, when available:

```text
GENCC_MOI
OMIM_INHERITANCE
CLINGEN_HI
CLINGEN_TS
ALLELE_GENOTYPE
SV_GENE_RELATIONSHIP
SVTYPE
```

The normalized fields include:

```text
GENE_MOI_SET
GENE_INHERITANCE_CLASS
SV_EFFECT_CLASS
INHERITANCE_MECHANISM_CLASS
INHERITANCE_MECHANISM_DETAIL
```

Examples of the intended interpretation:

- A DEL or direct disruptive breakpoint in a gene with **ClinGen HI=3** receives
  strong loss-of-function/dosage support.
- HI=2 is treated as emerging/supportive, not equivalent to HI=3.
- An AD disease model can support a direct loss/disruption candidate, but AD
  inheritance by itself does **not** prove that every possible SV mechanism is
  pathogenic. A heterozygous direct LoF/disruptive SV is treated as the clearest
  generic AD-compatible configuration. Unknown genotype is retained as an
  unresolved AD-compatible review state; a homozygous-alt SV in an AD gene is
  flagged for dedicated review rather than automatically promoted.
- Genes with both AD and AR disease models are assigned a
  `MIXED_AD_AR_DISEASE_MODEL_REVIEW` state unless disease-specific evidence
  resolves which model applies. The pipeline does not simply choose the more
  favorable inheritance mode.
- A whole-gene DUP is treated primarily as copy gain and is promoted by
  established/emerging **triplosensitivity (TS)** evidence, not simply because
  the gene has a dominant disease.
- A partial/intragenic DUP or a duplication breakpoint within a transcript is
  treated as a possible disruptive/LoF mechanism rather than automatically as
  pure copy gain. Orientation and breakpoint structure still require review.
- A direct heterozygous SV in an **AR gene** remains
  `DIRECT_EFFECT_RECESSIVE_SECOND_ALLELE_REQUIRED` unless a biallelic state or
  suitable second allele is demonstrated.
- X-linked events remain a sex/ploidy-specific review category.
- mtDNA genes remain a specialized review class because heteroplasmy and
  mtDNA-specific interpretation are not equivalent to nuclear AD/AR models.

ClinGen score 30 is treated as an AR association, not as evidence of
haploinsufficiency. Score 40 is retained as dosage-sensitivity-unlikely context.

## 4. Recessive SNV + SV pairing

The existing phased small-variant branch is used to look for a possible second
allele in AR genes.

`build_snv_sv_candidates.py` reports:

```text
AR_TRANS_SNV_SV_CANDIDATE
AR_CIS_NOT_BIALLELIC_BY_PHASE
AR_SECOND_ALLELE_CANDIDATE_PHASE_UNRESOLVED
NOT_AR_PAIRING_MODEL
```

Only a same-gene small variant phased in the same phase set and on the opposite
haplotype can be promoted to:

```text
AR_TRANS_SECOND_ALLELE_SUPPORTED
```

This supports a biallelic candidate model. It does not establish that either
allele is pathogenic and does not replace segregation/orthogonal validation.

## 5. Technical evidence

Technical evidence is no longer represented only by caller count.

The shared ranker uses:

```text
CALLER_COUNT
CALLER_READ_SUPPORT
CALLER_EVIDENCE_FLAGS
CALLER_EVIDENCE_MATCH
```

and reports:

```text
EVENT_TECHNICAL_TIER_V2
EVENT_MAX_CALLER_READ_SUPPORT
EVENT_TECHNICAL_REVIEW
```

The tiers are:

```text
STRONG
MODERATE
SUPPORTING
REVIEW
```

No caller-specific pathogenic weight is currently assigned by SV type. Such
weights would require benchmarking rather than being chosen arbitrarily.

## 6. Population frequency

needLR and exact gnomAD-SV remain distinct sources, but ranking uses the
**maximum explicit AF** across them. Therefore a common exact gnomAD-SV match
cannot be hidden by a lower/absent needLR value.

Default thresholds:

```text
ranking_rare_af: 0.001
ranking_max_af: 0.01
```

The main fields are:

```text
EVENT_POPULATION_TIER_V2
EVENT_MAX_EXPLICIT_AF
EVENT_MIN_EXPLICIT_AF
EVENT_POPULATION_AF_SOURCES
EVENT_POPULATION_AF_CONFLICT
```

Interpretation:

```text
VERY_RARE_MAX_AF_LE_0.001
LOW_FREQUENCY_MAX_AF_0.001_TO_0.01
UNKNOWN
TOO_COMMON_MAX_AF_GT_0.01
```

No match is not treated as AF=0. Missing population evidence is neutral.

AnnotSV benign-region AFmax is retained as overlap context but is not mixed into
the exact/provisional event AF because an overlapping benign SV is not
necessarily the same allele.

## 7. Constraint

LOEUF is used only as a modest late tie-break. Constraint is not a substitute
for disease mechanism.

Fields may include:

```text
GNOMAD_PLI
GNOMAD_LOEUF
GNOMAD_LOEUF_BIN
```

ClinGen HI/TS has stronger mechanistic meaning for dosage events and is handled
in the inheritance/mechanism layer.

## 8. Cohort recurrence

Cohort recurrence is an internal study signal, not population frequency.

```text
COHORT_SV_ID
COHORT_SAMPLE_COUNT
COHORT_RECURRENCE
```

An event found in more study samples is placed slightly lower only as a **late
cautionary tie-break**, after the major evidence axes. It is never removed,
because recurrence can represent a technical artifact, a common polymorphism,
a founder allele or a genuinely recurrent disease mechanism.

## 9. Panel versus non-panel ranking

Panel membership does not add points.

The pipeline creates explicit within-group ranks such as:

```text
GENE_RANK_WITHIN_PANEL_STATUS
EVENT_RANK_WITHIN_PANEL_STATUS
FINAL_GENE_RANK_WITHIN_PANEL_STATUS
FINAL_EVENT_RANK_WITHIN_PANEL_STATUS
MITO_RANK_WITHIN_ENCODING_PANEL_STATUS
```

Thesis ranking plots are written separately for panel and non-panel genes.

## 10. Mitochondrial ranking

Nuclear-encoded mitochondrial genes and mtDNA-encoded genes are not forced into
one common ranking.

Nuclear genes use the shared inheritance/mechanism, technical and population
logic. MitoCarta membership covers the broader mitochondrial proteome and is
not restricted to OXPHOS.

mtDNA genes are reported separately, including:

```text
MTDNA_PROTEIN_CODING
MTDNA_RRNA
MTDNA_TRNA
```

Their functional classes include OXPHOS complexes and mitochondrial translation
where applicable. mtDNA candidates still require mtDNA-specific interpretation,
including heteroplasmy, which is outside a generic nuclear SV rank.

## 11. Old fields retained for compatibility

Some older fields such as `PHENOTYPE_SCORE`,
`INTEGRATED_DISCOVERY_SCORE` and `ALLELE_RESEARCH_SCORE` remain in the
tables for provenance/backward compatibility.

They are **not** the primary final ordering model.

The current event ranking identifier is:

```text
sharedTieredRanking__geneTier_mechanismInheritance_technical_population_constraint__v6
```

The final patient-aware ranking may further incorporate patient HPO,
phased recessive pairing and cohort recurrence.

## 12. Ranking benchmark / sensitivity check

The repository includes:

```text
scripts/benchmark_candidate_ranking.py
```

It accepts a user-defined positive-control TSV containing `GENE` and,
optionally, `SV_ID`. It reports global and within-panel rank recovery,
top-k recall and rank percentile. This is intended for published/known positive
controls or deliberately planted test events; it does **not** fit weights.

Example:

```bash
python scripts/benchmark_candidate_ranking.py \
  --ranking SAMPLE_sv_gene_candidates.ranked.tsv \
  --truth known_positive_controls.tsv \
  --detail-output ranking_benchmark.detail.tsv \
  --summary-output ranking_benchmark.summary.tsv
```

This supports sensitivity analysis without training a model on the small thesis
cohort.

## 13. Interpretation limits

No field in these tables is a validated probability of pathogenicity.


The ranking does not replace:

- ACMG/AMP/ClinGen variant interpretation;
- segregation;
- trio analysis;
- breakpoint validation;
- functional experiments;
- mtDNA heteroplasmy analysis;
- clinical correlation.

The purpose is to order candidates transparently while preserving the reason a
candidate moved up or down.
