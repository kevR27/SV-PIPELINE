# LRS and SRS analysis parity

The two workflows follow the same biological reasoning, but they do not force
long-read and short-read data through identical software. A method is shared
when the same biological question can be answered reliably with both data
types. A method remains platform-specific when it depends on read length,
signal-level information, or a caller model trained for one technology.

## Shared biological analysis

| Biological question | LRS | SRS | Parity decision |
|---|---|---|---|
| Is the sequencing/alignment usable? | mosdepth, Dorado and read-length QC | mosdepth and samtools QC | Same question, platform-appropriate metrics |
| Which SVs are present genome-wide? | Sniffles2, cuteSV and DELLY, merged with Jasmine | Manta and DELLY, merged with SURVIVOR; CNVpytor adds read-depth CNVs | Same discovery principle; different callers |
| Is caller agreement visible without deleting single-caller events? | Complete Jasmine VCF plus multi-caller companion | Complete integrated SRS VCF plus multi-caller companion | Shared |
| Are small variants retained genome-wide? | Clair3 | DeepVariant | Same question; platform-specific caller |
| Are variants phased where the data allow it? | WhatsHap and LongPhase | WhatsHap | Shared core; LongPhase is LRS-specific |
| Are SVs annotated genome-wide? | AnnotSV and VEP | AnnotSV and VEP | Shared |
| Are panel and non-panel genes both retained? | Genome-wide table plus panel-only view | Genome-wide table plus panel-only view | Shared |
| Is gene relevance assessed with the same HON context? | Monarch/HPO and shared ranking scripts | Monarch/HPO and shared ranking scripts | Shared |
| Is allele-level evidence kept separate from pathogenicity? | Shared allele-assessment layer | Shared allele-assessment layer | Shared |
| Are nuclear mitochondrial genes and pathways identified? | MitoCarta and mitochondrial-gene ranking | MitoCarta and the same mitochondrial-gene ranking | Shared primary focus |
| Can top candidates be inspected in the reads? | Samplot long-read mode | Samplot short-read mode | Shared purpose; platform-specific display |

## Technology-specific analyses

| Analysis | Why it is not identical |
|---|---|
| LRS modkit methylation | Nanopore signal-level modification information is not present in ordinary Illumina reads. It is not expected in SRS. |
| LRS Straglr and TLDR | These methods use long-read span and sequence context. SRS uses ExpansionHunter and optional MELT for the corresponding variant classes. |
| LRS needLR | Its population backend and query preparation are based on long-read Sniffles calls. SRS must not reuse it. |
| SRS GRIDSS | Provides assembly/breakend support suited to paired short reads; it is not an extra requirement for the LRS callset. |
| SRS CNVpytor | Adds short-read depth and B-allele-frequency evidence. LRS uses its own coverage/context analyses instead. |
| SRS Mutserve2 | Provides the secondary short-read mtDNA SNV/heteroplasmy branch. It is not applied to Nanopore reads as if the error models were interchangeable. |

## Mitochondrial analysis hierarchy

The primary study question is whether a structural or small variant affects a
**nuclear gene that encodes a mitochondrial protein or participates in a
mitochondrial/bioenergetic pathway**. MitoCarta membership, pathway,
subcompartment, HON relevance, SV mechanism, technical evidence and population
evidence remain separate columns so the reason for prioritization is visible.

The mitochondrial genome is retained as a **secondary analysis**. In SRS,
Mutserve2 reports mtDNA SNV/heteroplasmy candidates separately. In LRS, chrM
small variants already present in the Clair3 VCF are copied to a separate
secondary review table, while chrM SVs can appear in the mitochondrial SV
ranking. The LRS Clair3 output is not labelled as validated heteroplasmy and is
not presented as a complete mtDNA diagnostic workflow.

## Shared source of truth

Neither workflow restricts discovery to the optic-neuropathy panel. The master
callset is genome-wide. Panel membership, nuclear mitochondrial biology and HPO
relevance are interpretation layers applied afterward. This preserves the
ability to investigate non-panel genes without presenting them as established
disease genes.
