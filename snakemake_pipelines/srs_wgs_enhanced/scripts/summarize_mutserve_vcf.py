#!/usr/bin/env python3
"""Create a compact mtDNA SNV/heteroplasmy table from Mutserve2 VCF output."""
import argparse,csv,gzip
from pathlib import Path

def openx(p):return gzip.open(p,"rt") if str(p).endswith(".gz") else open(p)
def main():
 p=argparse.ArgumentParser();p.add_argument("--vcf",required=True);p.add_argument("--output",required=True);a=p.parse_args()
 rows=[]
 with openx(a.vcf) as h:
  for line in h:
   if line.startswith("#"):continue
   f=line.rstrip().split("\t")
   if len(f)<8:continue
   info=dict(x.split("=",1) if "=" in x else (x,"True") for x in f[7].split(";") if x)
   fmt=dict(zip(f[8].split(":"),f[9].split(":"))) if len(f)>9 else {}
   rows.append({"CHROM":f[0],"POS":f[1],"ID":f[2],"REF":f[3],"ALT":f[4],"FILTER":f[6],
                "GT":fmt.get("GT","."),"HETEROPLASMY_AF":fmt.get("AF",info.get("AF",".")),
                "DP":fmt.get("DP",info.get("DP",".")),
                "NOTE":"Mutserve2 VCF is primarily mtDNA SNV/heteroplasmy evidence; its documentation notes indels are not included in VCF output."})
 out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
 cols=["CHROM","POS","ID","REF","ALT","FILTER","GT","HETEROPLASMY_AF","DP","NOTE"]
 with open(out,"w",newline="") as h:
  w=csv.DictWriter(h,fieldnames=cols,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
 print(f"[OK] mtDNA_records={len(rows)} output={out}")
if __name__=="__main__":main()
