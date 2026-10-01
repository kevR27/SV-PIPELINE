# HPO reference logic for hereditary optic neuropathy

## Active phenotype seed file

The active seed file is:

`reference/hon_hpo_seed_terms.tsv`

It contains only core hereditary optic neuropathy/optic-nerve phenotype anchors.
These terms are used to identify phenotype proximity in the Monarch knowledge
graph. They are not intended to enumerate every manifestation of every
syndromic optic neuropathy.

## Earlier manual term list

The earlier `reference/optic_neuropathy_hpo_terms.tsv` file has been removed
from the active repository. It mixed core ocular features with selected
non-ocular features without linking those features to the specific disease/gene
in which they occur. That structure is not appropriate for isolated-versus-
syndromic interpretation.

The only active manually curated seed file is
`reference/hon_hpo_seed_terms.tsv`. Broader syndromic phenotypes are obtained
from the generated gene-disease-HPO reference rather than from a second global
manual list.

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
