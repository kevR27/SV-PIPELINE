#!/usr/bin/env python3
"""Report coding/splice small variants in genes that also contain candidate SRS SVs.

This is an SRS-specific secondary review layer. It never infers cis/trans
between a small variant and an SV because the SV itself is not long-range
phased by the short-read Whatshap branch.
"""
from __future__ import annotations
import argparse,gzip,re
from pathlib import Path
import pandas as pd

KEEP_IMPACT={"HIGH","MODERATE"}
KEEP_CONSEQUENCE={"splice_acceptor_variant","splice_donor_variant"}

def open_text(path):
    return gzip.open(path,"rt") if str(path).endswith(".gz") else open(path,encoding="utf-8")

def read_phase(path):
    by_key={}; by_pos={}
    with open_text(path) as fh:
        sample_i=None
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                f=line.rstrip().split("\t"); sample_i=9 if len(f)>9 else None; continue
            if line.startswith("#") or sample_i is None: continue
            f=line.rstrip().split("\t")
            chrom=f[0].removeprefix("chr"); pos=int(f[1]); alts=f[4].split(",")
            fmt=f[8].split(":"); vals=f[sample_i].split(":"); d=dict(zip(fmt,vals))
            gt=d.get("GT","."); ps=d.get("PS",".")
            phased="YES" if "|" in gt and ps not in {"","."} else "NO"
            for alt in alts: by_key[(chrom,pos,alt)]=(gt,ps,phased)
            if len(alts)==1: by_pos[(chrom,pos)]=(gt,ps,phased)
    return by_key,by_pos

def read_vep(path):
    header=None; rows=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#Uploaded_variation"):
                header=line[1:].rstrip().split("\t"); continue
            if line.startswith("#") or header is None: continue
            rows.append(dict(zip(header,line.rstrip().split("\t"))))
    return rows

def extra_map(v):
    out={}
    for x in str(v or "").split(";"):
        if "=" in x:
            k,val=x.split("=",1); out[k]=val
    return out

def parse_uploaded(v):
    m=re.match(r"^(?:chr)?([^_]+)_(\d+)_([^/]+)/(.+)$",str(v))
    if not m: return None
    c,p,_r,a=m.groups(); return c,int(p),a

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--vep",required=True)
    ap.add_argument("--phased-vcf",required=True)
    ap.add_argument("--sv-candidates",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()

    sv=pd.read_csv(a.sv_candidates,sep="\t",dtype=str,low_memory=False)
    gene_col=next((x for x in ["GENE","Gene","ANNotsv_Gene","GENES"] if x in sv.columns),None)
    if gene_col is None: raise ValueError("SV candidate table has no gene column")
    sv["_GENE_NORMALIZED"]=sv[gene_col].fillna(".").astype(str)
    genes=set(g for value in sv["_GENE_NORMALIZED"] for g in re.split(r"[;,|/]",value) if g not in {"","."})

    phase,phase_pos=read_phase(a.phased_vcf)
    vep=read_vep(a.vep)
    small=[]
    for r in vep:
        ex=extra_map(r.get("Extra",""))
        gene=str(r.get("SYMBOL") or ex.get("SYMBOL","."))
        if gene not in genes: continue
        cons=str(r.get("Consequence","."))
        impact=str(r.get("IMPACT") or ex.get("IMPACT",".")).upper()
        if impact not in KEEP_IMPACT and not (set(cons.split(",")) & KEEP_CONSEQUENCE):
            continue
        u=parse_uploaded(r.get("Uploaded_variation","."))
        if not u: continue
        chrom,pos,alt=u
        gt,ps,phased=phase.get((chrom,pos,alt),phase_pos.get((chrom,pos),(".",".","NO")))
        small.append({
            "GENE":gene,"SMALL_VARIANT":r.get("Uploaded_variation","."),
            "CHROM":chrom,"POS":str(pos),"ALT":alt,
            "CONSEQUENCE":cons,"IMPACT":impact,
            "EXISTING_VARIATION":r.get("Existing_variation","."),
            "SMALL_GT":gt,"SMALL_PS":ps,"SMALL_PHASED":phased
        })

    rows=[]
    for sm in small:
        matches=sv[sv["_GENE_NORMALIZED"].apply(lambda x: sm["GENE"] in re.split(r"[;,|/]",str(x)))]
        for _,vr in matches.iterrows():
            rows.append({
                **sm,
                "SV_ID":vr.get("SV_ID","."),
                "SVTYPE":vr.get("SVTYPE","."),
                "SV_CHROM":vr.get("CHROM","."),
                "SV_START":vr.get("START","."),
                "SV_END":vr.get("END","."),
                "SV_GENE_RELATIONSHIP":vr.get("SV_GENE_RELATIONSHIP",vr.get("GENE_RELATIONSHIP",".")),
                "SV_GNOMAD_EXACT_MATCH":vr.get("GNOMAD_SV_EXACT_MATCH","."),
                "SV_GNOMAD_AF":vr.get("GNOMAD_SV_AF","."),
                "SV_CALLER_COUNT":vr.get("CALLER_COUNT","."),
                "SV_GRIDSS_SUPPORT":vr.get("GRIDSS_BREAKPOINT_MATCH","."),
                "SV_CNVPYTOR_SUPPORT":vr.get("CNVPYTOR_RD_MATCH","."),
                "PHASE_RELATION":"NOT_ESTABLISHED_FOR_SV_IN_SRS",
                "INTERPRETATION":"Same-gene small variant and SV; short-read small-variant phase does not establish cis/trans relative to the SV."
            })

    cols=["GENE","SMALL_VARIANT","CHROM","POS","ALT","CONSEQUENCE","IMPACT",
          "EXISTING_VARIATION","SMALL_GT","SMALL_PS","SMALL_PHASED","SV_ID",
          "SVTYPE","SV_CHROM","SV_START","SV_END","SV_GENE_RELATIONSHIP",
          "SV_GNOMAD_EXACT_MATCH","SV_GNOMAD_AF","SV_CALLER_COUNT",
          "SV_GRIDSS_SUPPORT","SV_CNVPYTOR_SUPPORT","PHASE_RELATION","INTERPRETATION"]
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows,columns=cols).to_csv(out,sep="\t",index=False)
    print(f"[OK] srs_snv_sv_pairs={len(rows)} output={out}")
if __name__=="__main__":
    main()
