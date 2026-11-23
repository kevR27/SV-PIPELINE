# HPO reference logic for hereditary optic neuropathy

## Active phenotype seed file

The active seed file is:

`reference/hon_hpo_seed_terms.tsv`

Its path is explicit in both configs under `hon_hpo_seed_terms`. The custom ranking and LRS HPO-reference generation use the same setting. The older `discovery_hpo_file` key remains an alias; conflicting
paths are rejected. `hpo_terms` accepts IDs, not a file path.

The current seed file has **6 rows**, all `CORE_OCULAR_HON`.
The older six generic `MITO_SYNDROMIC_CONTEXT` rows have been removed.

These terms are shared **ontology anchors**, not a statement that all six
findings have been observed in every patient. Individual patients' additional
systemic or neurological HPO terms must be entered explicitly into the
private `patient_profiles_tsv` table if clinically documented.

The supported workflow for patient-specific Genomiser analysis, age-onset
ranges and independently evaluated diagnostic hypotheses is documented in
`docs/GENOMISER_PATIENT_WORKFLOW.md`. Both the custom final ranking and the
Genomiser genome preset consume the *same present HPO terms* for each patient;
an absent term or suspected clinical diagnosis is not treated as a confirmed
positive HPO feature.

## Earlier manual term list

The earlier `reference/optic_neuropathy_hpo_terms.tsv` file has been removed
from the active repository. It mixed core ocular features with selected
non-ocular features without linking those features to the specific disease/gene
in which they occur. That structure is not appropriate for isolated-versus-
syndromic interpretation.

The only active manually curated seed file is
`reference/hon_hpo_seed_terms.tsv`. Its context rows have an explicit role
instead of being treated as shared patient findings. Broader disease-specific
syndromic phenotypes are obtained from the generated gene-disease-HPO reference
rather than from a second global manual list.

## Isolated and syndromic HON

Syndromic manifestations are generated from gene-disease-HPO associations in
the Monarch knowledge graph rather than inserted into one global HON term
list. The generated cohort reference is:

`cohort_analysis/reference/hon_hpo_reference.tsv`

It records HPO term, HPO name, HON relevance, disease, gene, MitoCarta pathway,
publication/knowledge-source provenance and descriptive phenotype breadth.

The pipeline does **not** automatically label a disease as isolated or
syndromic solely from HPO breadth. The current `ISOLATED_SYNDROMIC` field
therefore remains explicitly uncurated until a disease-level curated source or
manual literature curation is supplied.

## Relation to structural variants

HPO evidence is gene/disease-level phenotype evidence. It can prioritize a gene
overlapped or disrupted by an SV even when much of the historical literature
describes SNVs. However, that does not imply that every SV in the gene has the
same functional consequence as a reported SNV.

The phenotype layer and the SV mechanism layer must therefore remain separate:

`phenotype compatibility` + `SV-gene mechanism` + `population evidence`
+ `technical support` + `segregation/experimental validation when available`.

For example, an established missense gain-of-function disease mechanism does
not automatically support a whole-gene deletion. Conversely, a known
loss-of-function disease mechanism may be mechanistically compatible with a
gene-disrupting deletion or breakpoint, but still requires variant-specific
evidence.

## Patient-specific HPO

Patient-specific semantic similarity is optional until clinical HPO terms are
available. If no patient HPO file is configured, the pipeline reports
`PATIENT_HPO_NOT_AVAILABLE` and does not fabricate a patient-specific
similarity score.
