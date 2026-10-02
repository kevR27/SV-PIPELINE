#!/usr/bin/env python3
"""Generate Samplot images for prioritized Illumina SV candidates.

Unlike the LRS helper, this script uses Samplot in ordinary short-read mode and
never passes --long_read.
"""
from __future__ import annotations
import argparse,re,subprocess
from pathlib import Path
import pandas as pd

def first_col(df,names):
    low={str(c).lower():c for c in df.columns}
    for n in names:
        if n in df.columns: return n
        if n.lower() in low: return low[n.lower()]
    return None

def num(v):
    try: return int(float(v))
    except Exception: return None

def safe(v,n=55):
    x=re.sub(r"[^A-Za-z0-9_.-]+","_",str(v)).strip("_") or "candidate"
    return x[:n]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--bam",required=True)
    ap.add_argument("--reference",required=True)
    ap.add_argument("--sample",required=True)
    ap.add_argument("--out-dir",required=True)
    ap.add_argument("--top-n",type=int,default=20)
    ap.add_argument("--min-mapq",type=int,default=20)
    ap.add_argument("--window",type=int,default=5000)
    ap.add_argument("--zoom",type=int,default=20000)
    ap.add_argument("--large-sv-threshold",type=int,default=1000000)
    ap.add_argument("--gene-annotation")
    ap.add_argument("--manifest",required=True)
    a=ap.parse_args()

    df=pd.read_csv(a.input,sep="\t",dtype=str,low_memory=False)
    idc=first_col(df,["SV_ID","ID"]); cc=first_col(df,["CHROM"]); sc=first_col(df,["START","POS"])
    ec=first_col(df,["END"]); tc=first_col(df,["SVTYPE","SV_type"]); gc=first_col(df,["GENE","GENES","ANNotsv_Gene"])
    scorec=first_col(df,["EVENT_GENE_RELEVANCE_SCORE","INTEGRATED_DISCOVERY_SCORE","PHENOTYPE_SCORE","ALLELE_RESEARCH_SCORE"])
    callc=first_col(df,["CALLER_COUNT","SUPP"])
    if None in (idc,cc,sc,tc): raise ValueError("Input requires SV_ID, CHROM, START and SVTYPE")

    work=df.copy()
    work["_score"]=pd.to_numeric(work[scorec],errors="coerce").fillna(0) if scorec else 0
    work["_callers"]=pd.to_numeric(work[callc],errors="coerce").fillna(0) if callc else 0
    work["_gridss"]=work.get("GRIDSS_BREAKPOINT_MATCH",pd.Series("NO",index=work.index)).eq("YES").astype(int)
    work["_cnv"]=work.get("CNVPYTOR_RD_MATCH",pd.Series("NO",index=work.index)).eq("YES").astype(int)
    work=work.sort_values(["_score","_gridss","_cnv","_callers"],ascending=[False,False,False,False])
    work=work.drop_duplicates(idc).head(a.top_n)

    outdir=Path(a.out_dir); outdir.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for rank,(_,r) in enumerate(work.iterrows(),1):
        chrom=str(r[cc]); start=num(r[sc]); end=num(r[ec]) if ec else None
        if start is None: continue
        if end is None or end<=start: end=start+1
        svtype=str(r[tc]).upper()
        gene=str(r[gc]) if gc else "."
        sid=str(r[idc])
        out=outdir/f"{rank:02d}_{safe(svtype,12)}_{safe(chrom,12)}_{start}_{end}_{safe(sid)}_{safe(gene,30)}.png"
        cmd=["samplot","plot","-n",a.sample,"-b",a.bam,"-r",a.reference,"-o",str(out),
             "-c",chrom,"-s",str(start),"-e",str(end),"-q",str(a.min_mapq)]
        if svtype in {"DEL","DUP","INV"}: cmd += ["-t",svtype]
        if abs(end-start)>=a.large_sv_threshold: cmd += ["--zoom",str(a.zoom)]
        else: cmd += ["--window",str(a.window)]
        if a.gene_annotation:
            cmd += ["-A",a.gene_annotation,"--annotation_filenames","Genes","--annotation_fontsize","7"]
        p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        status="OK" if p.returncode==0 and out.exists() and out.stat().st_size>0 else "FAILED"
        manifest.append({
            "RANK":rank,"SV_ID":sid,"SVTYPE":svtype,"CHROM":chrom,"START":start,"END":end,
            "GENE":gene,"GRIDSS_SUPPORT":r.get("GRIDSS_BREAKPOINT_MATCH","."),
            "CNVPYTOR_SUPPORT":r.get("CNVPYTOR_RD_MATCH","."),
            "IMAGE":str(out) if status=="OK" else ".","STATUS":status,
            "MESSAGE":"." if status=="OK" else (p.stdout or "").strip().split("\n")[-1]
        })
    m=Path(a.manifest); m.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(manifest).to_csv(m,sep="\t",index=False)
    print(f"[OK] samplot_requested={len(work)} generated={sum(x['STATUS']=='OK' for x in manifest)}")
if __name__=="__main__":
    main()
