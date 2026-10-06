# Candidate-ranking interpretation

This pipeline uses **research prioritization**, not a single pathogenicity score.

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

## 1. Phenotype and gene relevance

The old anchor-count score (`min(10, anchor_count * 5)`) is no longer the
primary gene-relevance model.

Generic hereditary optic-neuropathy relevance is calculated with
**Resnik best-match-average semantic similarity** against the broader human
gene-HPO background. The reference is now split into two roles:

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
GENE_RELEVANCE_TIER
GENE_RELEVANCE_DISPLAY_SCORE
GENE_RELEVANCE_PHENOTYPE_SCOPE
```

The display score is balanced between phenotype and curated disease evidence
and is retained for interpretation/tie-breaking. The primary ordering uses the
broad relevance tier:

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
patient's Resnik BMA result:

```text
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
- A DUP is promoted by established/emerging **triplosensitivity (TS)** evidence,
  not simply because the gene has a dominant disease.
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

## 13. Interpretation limits

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
