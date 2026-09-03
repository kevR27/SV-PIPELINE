#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PKG="$ROOT/SV_PIPELINE_final_fix"

backup() { local f="$1"; [[ -f "$f" ]] && cp -n "$f" "$f.before_sv_pipeline_fix"; }

# Copy shared scripts/resources.
mkdir -p "$ROOT/scripts" "$ROOT/reference" "$ROOT/envs"
cp "$PKG/scripts/"*.py "$ROOT/scripts/"
cp "$PKG/reference/optic_neuropathy_hpo_terms.tsv" "$ROOT/reference/"
cp "$PKG/envs/monarch.yaml" "$ROOT/envs/"
cp "$PKG/envs/stranger.yaml" "$ROOT/envs/"
cp "$PKG/envs/tldr.yaml" "$ROOT/envs/"

# ----- SRS config normalization -----
SRS_CFG="$ROOT/snakemake_pipelines/srs/config_srs.yaml"
backup "$SRS_CFG"
python - "$SRS_CFG" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text()
repls={
    'bed: "/PANEL_OA/optic_neuropathy_panel.bed"':'bed: "PANEL_OA/optic_neuropathy_panel.bed"',
    'candidate_genes_list: "/PANEL_OA/optic_neuropathy_genes.txt"':'candidate_genes_list: "PANEL_OA/optic_neuropathy_genes.txt"',
    'exclude_bed: "/reference/delly_hg38_exclude.bed"':'exclude_bed: "reference/delly_hg38_exclude.bed"',
    'expansionhunter_catalog: "/reference/variant_catalog_hg38_repeats.json"':'expansionhunter_catalog: "reference/variant_catalog_hg38_repeats.json"',
    'discovery_hpo_file: "/reference/optic_neuropathy_hpo_terms.tsv"':'discovery_hpo_file: "reference/optic_neuropathy_hpo_terms.tsv"',
}
for a,b in repls.items(): s=s.replace(a,b)
p.write_text(s)
PY

# ----- LRS config normalization -----
LRS_CFG="$ROOT/snakemake_pipelines/lrs/config_lrs.yaml"
backup "$LRS_CFG"
python - "$LRS_CFG" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text()
s=s.replace('cpg: "chr1:1-100000"', 'cpg: null')
s=s.replace('vntr_bed: "/reference/human_GRCh38_no_alt_analysis_set.trf.bed"', 'vntr_bed: "reference/human_GRCh38_no_alt_analysis_set.trf.bed"')
s=s.replace('tldr_elements: "/reference/teref.ont.human.fa"', 'tldr_elements: "reference/teref.ont.human.fa"')
if 'monarch_base_url:' not in s: s += '\nmonarch_base_url: "https://api-v3.monarchinitiative.org/v3/api"\n'
if 'discovery_hpo_file:' not in s: s += 'discovery_hpo_file: "reference/optic_neuropathy_hpo_terms.tsv"\n'
if 'annotsv_benign_af:' not in s: s += 'annotsv_benign_af: 0.01\n'
if 'sniffles_mosaic:' not in s: s += 'sniffles_mosaic: false\n'
p.write_text(s)
PY

# ----- LRS Snakefile deterministic patch -----
LRS="$ROOT/snakemake_pipelines/lrs/Snakefile_LRS"
backup "$LRS"
python - "$LRS" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text()
# Config-derived parameters
anchor='NEEDLR_BACKEND = config.get("needlr_backend_files", None)\n'
if 'MONARCH_BASE_URL = config.get(' not in s:
    s=s.replace(anchor, anchor + '\nMONARCH_BASE_URL = config.get("monarch_base_url", "https://api-v3.monarchinitiative.org/v3/api")\nDISCOVERY_HPO_FILE = config.get("discovery_hpo_file", "reference/optic_neuropathy_hpo_terms.tsv")\nANNOTSV_BENIGN_AF = config.get("annotsv_benign_af", 0.01)\nSNIFFLES_MOSAIC = config.get("sniffles_mosaic", False)\nSNIFFLES_MOSAIC_FLAG = "--mosaic" if SNIFFLES_MOSAIC else ""\nCPG_REGION = config.get("cpg", None)\nCPG_FLAG = f"--region {CPG_REGION}" if CPG_REGION else ""\nSTRAGLR_SCRIPT = config.get("straglr_script", "straglr.py")\n')
# rule all: add gene discovery, and remove Stranger from default target
s=s.replace('        expand(PATH + "{sample}/sv/straglr/{sample}_straglr.annotated.vcf", sample=SAMPLES),\n', '')
marker='        expand(PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.panel_only.tsv", sample=SAMPLES),\n'
add='''        expand(PATH + "{sample}/gene_discovery/{sample}_SV_genes.txt", sample=SAMPLES),\n        expand(PATH + "{sample}/gene_discovery/{sample}_nonpanel_genes.txt", sample=SAMPLES),\n        expand(PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv", sample=SAMPLES),\n        expand(PATH + "{sample}/gene_discovery/{sample}_ranked_candidates.tsv", sample=SAMPLES),\n'''
if 'gene_discovery/{sample}_ranked_candidates.tsv' not in s:
    s=s.replace(marker, marker+add)
# Sniffles mosaic made configurable
s=s.replace('            --minsupport 3 --no-qc --mosaic --allow-overwrite --min-alignment-length 100 \\\n', '            --minsupport 3 --no-qc {SNIFFLES_MOSAIC_FLAG} --allow-overwrite --min-alignment-length 100 \\\n')
# insertion FASTA: skip symbolic ALT records
s=s.replace("        bcftools view -i 'INFO/SVTYPE=\"INS\"' {input.vcf} | \\\n        bcftools query -f '>%ID\\n%ALT\\n' - | \\\n        sed 's/,.*$//' > {output.fasta}", "        bcftools view -i 'INFO/SVTYPE=\"INS\"' {input.vcf} | \\\n        bcftools query -f '>%ID\\t%ALT\\n' - | \\\n        awk -F'\\t' 'substr($2,1,1)!="<" && $2!="." {print \">"$1"\\n"$2}' > {output.fasta}")
# remove invalid Stranger rule block entirely
start=s.find('rule sv_annotation_stranger:')
if start!=-1:
    end=s.find('#####################################################################\n# needLR', start)
    if end!=-1:
        s=s[:start] + s[end:]
# needLR copy robustness + index
s=s.replace('''        cp {params.outdir}/*/*_RESULTS.tsv {output.tsv}\n        cp {params.outdir}/*/*_RESULTS.vcf.gz {output.vcf}''', '''        results_tsv=$(find {params.outdir} -type f -name '*_RESULTS.tsv' -print -quit)\n        results_vcf=$(find {params.outdir} -type f -name '*_RESULTS.vcf.gz' -print -quit)\n        test -n "$results_tsv"\n        test -n "$results_vcf"\n        cp "$results_tsv" {output.tsv}\n        cp "$results_vcf" {output.vcf}\n        tabix -f -p vcf {output.vcf}''')
# needLR description
s=s.replace('# needLR — population-AF filter ONLY (ONT-native 1KGP-LRSC cohort).', '# needLR — ONT-native population/reference SV annotation and filtering.')
# AnnotSV benign AF
if '-benignAF {ANNOTSV_BENIGN_AF}' not in s:
    s=s.replace('            -candidateGenesFiltering 0 \\\n            {params.hpo_flag}', '            -candidateGenesFiltering 0 \\\n            -benignAF {ANNOTSV_BENIGN_AF} \\\n            {params.hpo_flag}')
# Replace panel-only rule
start=s.find('rule annotsv_annotation_panel_only:')
end=s.find('##################################################################\n# VEP', start)
if start!=-1 and end!=-1:
    panel='''rule annotsv_annotation_panel_only:\n    input:\n        vcf = rules.needlr_annotation.output.vcf,\n        annotations = rules.annotsv_install_annotations.output.marker,\n        genes = GENE_LIST\n    output:\n        tsv = PATH + "{sample}/sv/annotsv/{sample}_merged_SV.annotsv.panel_only.tsv"\n    params:\n        outdir = PATH + "{sample}/sv/annotsv"\n    conda:\n        CONDAENV + "annotsv.yaml"\n    shell:\n        """\n        mkdir -p {params.outdir}\n        AnnotSV -SVinputFile {input.vcf} \\\n            -annotationsDir {ANNOTSV_ANNOTATIONS_DIR} \\\n            -outputDir {params.outdir} \\\n            -outputFile {wildcards.sample}_merged_SV.annotsv.panel_only.tsv \\\n            -genomeBuild GRCh38 \\\n            -candidateGenesFile {input.genes} \\\n            -candidateGenesFiltering 1 \\\n            -benignAF {ANNOTSV_BENIGN_AF}\n        """\n#####################################\n# Genome-wide human gene discovery\n#####################################\nrule extract_sv_genes:\n    input:\n        annotsv=rules.annotsv_annotation.output.tsv\n    output:\n        genes=PATH + "{sample}/gene_discovery/{sample}_SV_genes.txt"\n    conda:\n        CONDAENV + "monarch.yaml"\n    shell:\n        """\n        mkdir -p {PATH}/{wildcards.sample}/gene_discovery\n        python scripts/extract_annotsv_genes.py --annotsv {input.annotsv} --output {output.genes}\n        """\n\nrule extract_nonpanel_genes:\n    input:\n        genes=rules.extract_sv_genes.output.genes,\n        panel=GENE_LIST\n    output:\n        genes=PATH + "{sample}/gene_discovery/{sample}_nonpanel_genes.txt"\n    conda:\n        CONDAENV + "monarch.yaml"\n    shell:\n        """\n        python scripts/extract_nonpanel_genes.py --genes {input.genes} --panel {input.panel} --output {output.genes}\n        """\n\nrule human_gene_phenotypes:\n    input:\n        genes=rules.extract_nonpanel_genes.output.genes,\n        hpo=DISCOVERY_HPO_FILE\n    output:\n        tsv=PATH + "{sample}/gene_discovery/{sample}_human_gene_phenotypes.tsv"\n    params:\n        api=MONARCH_BASE_URL\n    conda:\n        CONDAENV + "monarch.yaml"\n    shell:\n        """\n        python scripts/monarch_human_gene_phenotypes.py \\\n            --genes {input.genes} --hpo-anchors {input.hpo} \\\n            --output {output.tsv} --api-base {params.api}\n        """\n\nrule rank_genomewide_candidates:\n    input:\n        annotsv=rules.annotsv_annotation.output.tsv,\n        genes=rules.extract_nonpanel_genes.output.genes,\n        phenotypes=rules.human_gene_phenotypes.output.tsv,\n        panel=GENE_LIST\n    output:\n        tsv=PATH + "{sample}/gene_discovery/{sample}_ranked_candidates.tsv"\n    conda:\n        CONDAENV + "monarch.yaml"\n    shell:\n        """\n        python scripts/rank_sv_gene_candidates.py \\\n            --annotsv {input.annotsv} --genes {input.genes} \\\n            --phenotypes {input.phenotypes} --panel {input.panel} \\\n            --output {output.tsv}\n        """\n'''
    s=s[:start]+panel+s[end:]
# modkit optional region
s=s.replace('{MODKITENV} extract calls -t {THREADS} --reference {input.ref} --cpg --region {CPG} {input.bam} {output.tsv}', '{MODKITENV} extract calls -t {THREADS} --reference {input.ref} --cpg {CPG_FLAG} {input.bam} {output.tsv}')
# LongPhase: merged SVs + ONT
s=s.replace('        sv_vcf = rules.sniffles2_sv.output.vcf,', '        sv_vcf = rules.jasmine_merge.output.vcf,')
s=s.replace('            longphase phase -s {input.vcf} --sv-file {input.sv_vcf} \\\n', '            longphase phase -s {input.vcf} --sv-file {input.sv_vcf} \\\n')
s=s.replace('            -b {input.bam} -r {input.ref} -t {THREADS} -o {params.prefix}', '            -b {input.bam} -r {input.ref} -t {THREADS} -o {params.prefix} --ont')
# add configurable Straglr script
s=s.replace('python /path/to/straglr/straglr.py', '{STRAGLR_SCRIPT}')
p.write_text(s)
PY

# ----- Add high-confidence SRS companion VCF without losing union callset -----
SRS="$ROOT/snakemake_pipelines/srs/Snakefile_SRS"
backup "$SRS"
python - "$SRS" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); s=p.read_text()
# Keep existing SURVIVOR union as source of truth; add high-confidence SUPP>=2 output.
if 'survivor_high_confidence' not in s:
    marker='#####################################\n# AnnotSV\n#####################################\n'
    block='''#####################################\n# High-confidence SV companion set\n#####################################\nrule survivor_high_confidence:\n    input:\n        vcf=rules.survivor_merge.output.vcf\n    output:\n        vcf=PATH + "{sample}/sv/merged/{sample}_merged_SV.SUPP2.vcf"\n    conda:\n        CONDAENV + "annotsv.yaml"\n    shell:\n        """\n        bcftools view -i 'INFO/SUPP>=2' -Ov {input.vcf} -o {output.vcf}\n        """\n'''
    s=s.replace(marker, block+marker, 1)
# Make discovery use existing genome-wide AnnotSV unchanged.
p.write_text(s)
PY

# Make a validation script.
cat > /mnt/data/SV_PIPELINE_final_fix/scripts/validate_setup.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys, yaml, re
root=Path(__file__).resolve().parents[1]
errors=[]
for p in [root/'snakemake_pipelines/srs/Snakefile_SRS', root/'snakemake_pipelines/lrs/Snakefile_LRS', root/'snakemake_pipelines/srs/config_srs.yaml', root/'snakemake_pipelines/lrs/config_lrs.yaml']:
    if not p.exists(): errors.append(f'Missing: {p}')
for p in sorted((root/'envs').glob('*.yaml')):
    try: yaml.safe_load(p.read_text())
    except Exception as e: errors.append(f'Invalid YAML: {p}: {e}')
for p in (root/'scripts').glob('*.py'):
    r=subprocess.run([sys.executable,'-m','py_compile',str(p)],capture_output=True,text=True)
    if r.returncode: errors.append(f'Python syntax error: {p}: {r.stderr.strip()}')
for snake in [root/'snakemake_pipelines/srs/Snakefile_SRS',root/'snakemake_pipelines/lrs/Snakefile_LRS']:
    if snake.exists() and '/path/to/' in snake.read_text(): errors.append(f'Placeholder path remains: {snake}')
if errors:
    print('FAIL')
    print('\n'.join('- '+e for e in errors)); sys.exit(1)
print('PASS: YAML parsing and Python syntax checks completed; Snakefiles contain no /path/to/ placeholders.')
print('Run Snakemake --lint and -n in the real analysis environment for full DAG validation.')
PY
chmod +x /mnt/data/SV_PIPELINE_final_fix/scripts/validate_setup.py /mnt/data/SV_PIPELINE_final_fix/patches/APPLY_CHANGES.sh

# Documentation
cat > /mnt/data/SV_PIPELINE_final_fix/README.md <<'EOF'
# Final SRS + LRS pipeline fixes

This package applies the remaining architecture/tooling fixes while preserving platform-specific SV calling.

## SRS stays unique
Manta + Delly -> SURVIVOR -> genome-wide AnnotSV.
The union callset remains the source of truth. A separate SUPP>=2 companion VCF is added for high-confidence review without discarding single-caller discovery calls.

## LRS stays unique
Sniffles2 + cuteSV -> Jasmine -> needLR -> AnnotSV.
The LRS fixes include:
- real panel-only AnnotSV run (`candidateGenesFiltering 1`);
- shared genome-wide gene-discovery branch;
- configurable Sniffles2 mosaic mode, default false;
- needLR output handling made robust;
- needLR documented correctly as population/reference annotation + filtering;
- removal of the invalid Straglr -> Stranger default path (Stranger is for ExpansionHunter/TRGT-style STR outputs);
- Straglr uses the configured executable rather than `/path/to/...`;
- LongPhase receives the Jasmine SV callset and `--ont`;
- modkit no longer runs against a hard-coded chr1 placeholder region;
- symbolic insertion alleles are excluded before RepeatMasker.

## Shared discovery
Both pipelines converge after AnnotSV:

SVs -> AnnotSV -> SV-associated genes -> remove optic-neuropathy panel -> human HGNC/HPO associations -> candidate ranking.

The HPO layer is human-only and represents gene-associated phenotype evidence, not a prediction of the patient's phenotype. Ranked candidates are not automatically classified as pathogenic.

## Apply
From the repository root, copy this directory to `SV_PIPELINE_final_fix`, then run:

`bash SV_PIPELINE_final_fix/patches/APPLY_CHANGES.sh`

The script creates `.before_sv_pipeline_fix` backups before modifying SRS/LRS config/Snakefiles.

## Validate
`python scripts/validate_setup.py`

Then:

`snakemake -s snakemake_pipelines/srs/Snakefile_SRS --lint`

`snakemake -s snakemake_pipelines/lrs/Snakefile_LRS --lint`

`snakemake -s snakemake_pipelines/srs/Snakefile_SRS -n --cores 1`

`snakemake -s snakemake_pipelines/lrs/Snakefile_LRS -n --cores 1`

External paths/databases and actual tool installations still need to exist on the analysis server.
