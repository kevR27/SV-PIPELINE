#!/usr/bin/env python3
"""Fail early on resources required by the isolated SRS-WGS workflow."""
from __future__ import annotations
import argparse,os
from pathlib import Path

def require(path,label):
    p=Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reference",required=True)
    ap.add_argument("--panel-bed",required=True)
    ap.add_argument("--gene-list",required=True)
    ap.add_argument("--gridss-enabled",choices=["true","false"],default="false")
    ap.add_argument("--mt-enabled",choices=["true","false"],default="false")
    ap.add_argument("--mt-contig",default="chrM")
    ap.add_argument("--bam",action="append",default=[])
    ap.add_argument("--bai",action="append",default=[])
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    ref=require(a.reference,"reference")
    require(str(ref)+".fai","reference FASTA index")
    require(a.panel_bed,"panel BED")
    require(a.gene_list,"candidate gene list")
    for b in a.bam: require(b,"BAM")
    for b in a.bai: require(b,"BAM index")
    if a.gridss_enabled=="true":
        missing=[str(ref)+s for s in (".amb",".ann",".bwt",".pac",".sa") if not Path(str(ref)+s).exists()]
        if missing:
            raise FileNotFoundError("GRIDSS default bwa backend requires reference bwa indexes: "+", ".join(missing))
    if a.mt_enabled=="true":
        contigs=set()
        with open(str(ref)+".fai",encoding="utf-8") as fh:
            for line in fh: contigs.add(line.split("\t",1)[0])
        if a.mt_contig not in contigs:
            raise ValueError(f"mtDNA contig {a.mt_contig!r} is absent from {ref}.fai")
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text("SRS_WGS_V2_PREFLIGHT_OK\n",encoding="utf-8")
if __name__=="__main__":
    main()
