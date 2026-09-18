# LRS pipeline logic

This document describes the active logic of `Snakefile_LRS_update`.

## Source-of-truth hierarchy

The per-patient Jasmine VCF is the master structural-variant (SV) callset. A downstream annotation tool is not allowed to remove a master SV simply because that tool cannot evaluate it.

```text
ONT BAM
├── mosdepth                         QC
├── Clair3 -> WhatsHap               complementary panel SNV/indel + phasing
├── Sniffles2 ─┐
├── cuteSV ────┼─ caller evidence QC -> Jasmine -> MASTER SV VCF
└── DELLY ─────┘                           ├── caller-support summary
                                          ├── AnnotSV -> panel/full views
                                          ├── VEP
                                          ├── gene + Monarch/HPO discovery
                                          └── integrated SV/gene table

ONT BAM -> dedicated Sniffles2 v2.6.2 query -> needLR -> population-frequency evidence
                                                       └-> matched back to MASTER SVs

ONT BAM -> Straglr -> repeat-expansion evidence
ONT BAM -> modkit extract -> per-read modification evidence
ONT BAM -> TLDR -> mobile-element insertion evidence
Clair3 + same-patient Sniffles SV VCF + BAM -> LongPhase
```

## Jasmine

The three genome-wide SV callers are normalized and transparently filtered before Jasmine. The input order is fixed as:

1. Sniffles2
2. cuteSV
3. DELLY

That order defines the interpretation of `SUPP_VEC`.

Jasmine is run per patient with `--allow_intrasample`. `--output_genotypes` is deliberately not used because these are three caller representations of the same biological sample, not three independent samples.

Jasmine first writes a plain raw VCF. The workflow then repairs missing metadata, sorts with `bcftools`, BGZF-compresses, and tabix-indexes the master VCF.

The complete Jasmine VCF remains the master callset. A separate high-confidence companion VCF is generated using the configured minimum caller count; it does not replace the master callset.

## needLR

needLR is a supplementary population-frequency branch, not a gatekeeper for AnnotSV.

A dedicated Sniffles2 v2.6.2 query is generated with the parameterization recommended by needLR as closely as possible. The cross-caller Jasmine consensus is not used as the query because the needLR backend was built from Sniffles2 calls.

The active workflow consumes `*_RESULTS.tsv` for population-frequency evidence. The needLR result VCF is retained as a secondary artifact but is not required to be tabix-indexed. This avoids failing an otherwise successful needLR analysis solely because its delivered VCF is not tabix-ready.

BNDs and SVs >=10 Mb remain in the Jasmine/AnnotSV master analysis. The integrated table marks them as not evaluable by needLR rather than treating missing needLR evidence as AF=0 or as evidence of rarity.

## AnnotSV and VEP

AnnotSV and VEP both consume the complete Jasmine master VCF directly. They do not consume a needLR-derived or needLR-filtered VCF.

The panel result is derived after genome-wide annotation. Panel membership is therefore an annotation/view, not a calling restriction.

## Integrated evidence table

`build_integrated_sv_gene_tsv.py` defines rows from the Jasmine master callset and combines:

- Jasmine caller provenance and `SUPP_VEC`
- filtered caller evidence from Sniffles2, cuteSV, and DELLY
- AnnotSV gene/clinical annotation
- needLR population-frequency evidence matched by SV type and genomic location
- panel membership
- Monarch/HPO phenotype evidence

needLR evidence must never be joined merely because two SVs affect the same gene.

## Legacy helpers

The repository still contains `prepare_needlr_vcf.py` and `fix_needlrvcf_for_annotsv.py` for historical/reference purposes. They are **not part of the active `Snakefile_LRS_update` workflow**. Do not rebuild a `Jasmine -> needLR -> repaired VCF -> AnnotSV` chain with them.

## Methylation note

The active workflow currently keeps `modkit extract calls`. A `modkit pileup`/bedMethyl branch is intentionally not enabled in this revision because current modkit 0.6.4 has open pileup correctness/stability issues on some valid inputs. Revisit bedMethyl generation after an upstream release addresses those issues and after validation on the project modBAMs.

## Recommended validation order

Before running the cohort, validate one known patient end-to-end:

1. Parse and filter each caller.
2. Confirm DELLY retains calls and is not systematically marked `LOW_SUPPORT`.
3. Produce the Jasmine master VCF and inspect `SUPP`, `SUPP_VEC`, and `IDLIST`.
4. Run the dedicated Sniffles2 needLR query and confirm `*_RESULTS.tsv` is produced.
5. Run AnnotSV directly on the Jasmine master VCF.
6. Build the integrated table and inspect `NEEDLR_STATUS` rather than interpreting missing needLR AF as zero.
7. Confirm a known positive SV survives calling -> evidence QC -> Jasmine -> AnnotSV -> integrated table.

Only after the single-patient validation should the same workflow be expanded across the cohort.
