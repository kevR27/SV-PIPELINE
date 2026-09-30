# Inversion interpretation

The pipeline does not treat every gene inside an inversion as disrupted.

An inversion is usually copy-number neutral. Its main direct effects are expected
at the breakpoints, although regulatory position effects can occur when a
breakpoint changes enhancer-promoter relationships or topological domains.

The final SV-gene table therefore uses the following hierarchy:

| SV_GENE_EFFECT | Interpretation |
| --- | --- |
| `INVERSION_TWO_BREAKPOINTS_IN_GENE` | Both inversion breakpoints intersect the same gene. High-priority structural review. |
| `INVERSION_BREAKPOINT_IN_EXON` | A breakpoint intersects an exon. Direct coding disruption is plausible. |
| `INVERSION_BREAKPOINT_IN_INTRON` | A breakpoint intersects the transcribed gene in an intron. Gene disruption/splicing effects are possible but not assumed. |
| `INVERSION_BREAKPOINT_IN_TRANSCRIPT` | The breakpoint intersects the transcript but the exon/intron context is unresolved. |
| `INVERSION_BREAKPOINT_NEAR_GENE` | Breakpoint lies within the configured proximity window. This is a regulatory/position-effect hypothesis, not direct disruption. |
| `GENE_INSIDE_INVERSION` | The gene lies inside the inverted segment but neither breakpoint intersects or approaches it. This is retained as interval context only. |
| `INVERSION_GENE_EFFECT_UNRESOLVED` | The available annotation is insufficient to resolve the relationship. |

The default proximity window is 10 kb. This is a research-prioritization window,
not a pathogenicity threshold.

For a candidate inversion, review should proceed in this order:

1. verify that the inversion itself is technically credible from caller/read evidence;
2. inspect both breakpoints in Samplot/IGV;
3. determine whether either breakpoint directly intersects a relevant transcript;
4. check for small deletion/duplication imbalance around the breakpoints;
5. if the breakpoints are noncoding, consider promoter/enhancer/TAD position effects;
6. review population/internal recurrence and, when available, inheritance;
7. use RNA or other functional evidence only when biologically appropriate.

Large inversions should therefore not inflate gene ranking simply because they
contain hundreds or thousands of genes.

Relevant examples and approaches:
- PhenoSV evaluates inversion breakpoints separately rather than assigning the
  full inverted interval as equivalent gene disruption:
  https://www.nature.com/articles/s41467-023-43651-y
- Review of disease mechanisms caused by structural disruption of 3D genome
  architecture:
  https://www.nature.com/articles/s41576-025-00862-x
- Long-read resolution of rare/pathogenic inversions:
  https://pubmed.ncbi.nlm.nih.gov/39486878/
