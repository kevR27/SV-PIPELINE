#!/usr/bin/env python3
"""Parse CNVpytor text calls from one or more bin sizes.

CNVpytor's documented -call output begins with:
  CNV type, region, size, normalized RD level, followed by statistical fields.

This parser deliberately preserves the raw row and does not reinterpret the
statistics beyond the documented first fields. The primary-bin calls are also
exported as a simple symbolic-SV VCF for AnnotSV.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

MISSING="."

def parse_region(value: str):
    m=re.match(r"([^:]+):(\d+)-(\d+)$", value.replace(",",""))
    if not m:
        raise ValueError(f"Unrecognized CNVpytor region: {value}")
    chrom,start,end=m.groups()
    return chrom,int(start),int(end)

def normalize_type(value: str):
    v=value.lower()
    if v.startswith("del"):
        return "DEL"
    if v.startswith("dup"):
        return "DUP"
    return v.upper()

def safe_float(value):
    try:
        return float(value)
    except Exception:
        return None

def read_calls(bin_size: int, path: str):
    rows=[]
    with open(path,"r",encoding="utf-8",errors="replace") as fh:
        for line_no,line in enumerate(fh,1):
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            fields=re.split(r"\s+",line)
            if len(fields)<4:
                continue
            try:
                chrom,start,end=parse_region(fields[1])
            except ValueError:
                continue
            svtype=normalize_type(fields[0])
            if svtype not in {"DEL","DUP"}:
                continue
            size=fields[2]
            rd_level=fields[3]
            event_id=f"CNVPYTOR_{bin_size}_{chrom}_{start}_{end}_{svtype}"
            rows.append({
                "CNVPYTOR_ID":event_id,
                "BIN_SIZE":str(bin_size),
                "CHROM":chrom,
                "START":str(start),
                "END":str(end),
                "SVTYPE":svtype,
                "SVLEN":str(-(end-start) if svtype=="DEL" else (end-start)),
                "REPORTED_SIZE":size,
                "RD_LEVEL":rd_level,
                "EVAL1":fields[4] if len(fields)>4 else MISSING,
                "EVAL2":fields[5] if len(fields)>5 else MISSING,
                "EVAL3":fields[6] if len(fields)>6 else MISSING,
                "EVAL4":fields[7] if len(fields)>7 else MISSING,
                "Q0":fields[8] if len(fields)>8 else MISSING,
                "PN":fields[9] if len(fields)>9 else MISSING,
                "DG":fields[10] if len(fields)>10 else MISSING,
                "RAW_CALL":line,
            })
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--call", action="append", required=True,
                    help="BIN=path; repeat for each CNVpytor bin")
    ap.add_argument("--primary-bin", type=int, required=True)
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--output-vcf", required=True)
    args=ap.parse_args()

    all_rows=[]
    for item in args.call:
        if "=" not in item:
            raise ValueError("--call must be BIN=path")
        b,p=item.split("=",1)
        all_rows.extend(read_calls(int(b),p))

    out=Path(args.output_tsv)
    out.parent.mkdir(parents=True,exist_ok=True)
    cols=["CNVPYTOR_ID","BIN_SIZE","CHROM","START","END","SVTYPE","SVLEN",
          "REPORTED_SIZE","RD_LEVEL","EVAL1","EVAL2","EVAL3","EVAL4","Q0","PN","DG","RAW_CALL"]
    with out.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(all_rows)

    primary=[r for r in all_rows if int(r["BIN_SIZE"])==args.primary_bin]
    vcf=Path(args.output_vcf)
    vcf.parent.mkdir(parents=True,exist_ok=True)
    contigs=[]
    seen=set()
    for r in primary:
        if r["CHROM"] not in seen:
            seen.add(r["CHROM"]); contigs.append(r["CHROM"])
    with vcf.open("w",encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=CNVpytor\n")
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n')
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">\n')
        fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        fh.write('##INFO=<ID=CNVPYTOR_BIN,Number=1,Type=Integer,Description="CNVpytor bin size">\n')
        fh.write('##INFO=<ID=RD_LEVEL,Number=1,Type=Float,Description="CNVpytor normalized read-depth level">\n')
        for c in contigs:
            fh.write(f"##contig=<ID={c}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in primary:
            info=f"END={r['END']};SVTYPE={r['SVTYPE']};SVLEN={r['SVLEN']};CNVPYTOR_BIN={r['BIN_SIZE']};RD_LEVEL={r['RD_LEVEL']}"
            fh.write(f"{r['CHROM']}\t{r['START']}\t{r['CNVPYTOR_ID']}\tN\t<{r['SVTYPE']}>\t.\tPASS\t{info}\n")

    print(f"[OK] CNVpytor rows={len(all_rows)} primary={len(primary)}")

if __name__=="__main__":
    main()
