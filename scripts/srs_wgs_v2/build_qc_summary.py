#!/usr/bin/env python3
"""Summarize stable fields from samtools flagstat/stats and mosdepth."""
from __future__ import annotations
import argparse,csv,re
from pathlib import Path

def parse_flagstat(path):
    out={}
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            m=re.match(r"(\d+) \+ (\d+) (.+)",line.strip())
            if not m: continue
            passed,failed,label=m.groups()
            key=re.sub(r"[^A-Za-z0-9]+","_",label.split("(")[0]).strip("_").lower()
            out[f"flagstat_{key}_qcpass"]=passed
            out[f"flagstat_{key}_qcfail"]=failed
    return out

def parse_stats(path):
    out={}
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if not line.startswith("SN\t"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<3: continue
            key=re.sub(r"[^A-Za-z0-9]+","_",f[1].rstrip(":")).strip("_").lower()
            out[f"samtools_{key}"]=f[2]
    return out

def parse_mosdepth(path):
    rows=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        rd=csv.DictReader(fh,delimiter="\t")
        rows=list(rd)
    out={}
    for r in rows:
        chrom=r.get("chrom","")
        if chrom=="total":
            for k,v in r.items():
                if k!="chrom": out[f"mosdepth_total_{k}"]=v
    if not out and rows:
        r=rows[-1]
        for k,v in r.items():
            if k!="chrom": out[f"mosdepth_lastrow_{k}"]=v
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample",required=True)
    ap.add_argument("--flagstat",required=True)
    ap.add_argument("--stats",required=True)
    ap.add_argument("--mosdepth-summary",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    row={"sample":a.sample}
    row.update(parse_flagstat(a.flagstat))
    row.update(parse_stats(a.stats))
    row.update(parse_mosdepth(a.mosdepth_summary))
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(row),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerow(row)
if __name__=="__main__":
    main()
