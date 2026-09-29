# MitoCarta3.0 reference

This pipeline does not vendor the MitoCarta workbook. Keep the original external
reference file so provenance remains explicit.

Download the official Human MitoCarta3.0 workbook from the Broad Institute:

```bash
mkdir -p reference/mitocarta

wget \
  https://personal.broadinstitute.org/scalvo/MitoCarta_Download/Human.MitoCarta3.0.xls \
  -O reference/mitocarta/Human.MitoCarta3.0.xls
```

The LRS post-processing workflow reads the workbook from `mitocarta_file` in
`config_lrs.yaml`.

The annotation layer adds MitoCarta membership, nuclear-versus-mtDNA encoding,
Maestro score, localization evidence, sub-mitochondrial compartment,
MitoPathways, top-level mitochondrial pathways, and a combined
`MITO_ON_ASSOCIATION_CLASS`.

`MITO_ON_ASSOCIATION_CLASS` is a descriptive co-annotation:
a gene is both MitoCarta-positive and linked to the optic-neuropathy phenotype
context already computed by the pipeline. It is not a correlation coefficient,
a pathogenicity probability, or proof that the SV causes the phenotype.

Source: MitoCarta3.0, Broad Institute / Rath et al., Nucleic Acids Research.
