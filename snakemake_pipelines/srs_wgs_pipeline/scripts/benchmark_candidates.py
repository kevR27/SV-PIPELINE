#!/usr/bin/env python3
"""Measure candidate-list recall and rank against a positive-control TSV."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd


def find(df, names):
    lower = {str(c).lower(): c for c in df.columns}
    return next((n if n in df.columns else lower[n.lower()] for n in names if n in df.columns or n.lower() in lower), None)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("--candidates", required=True); p.add_argument("--truth", required=True)
    p.add_argument("--output-tsv", required=True); p.add_argument("--output-json", required=True); a = p.parse_args()
    candidates = pd.read_csv(a.candidates, sep="\t", dtype=str); truth = pd.read_csv(a.truth, sep="\t", dtype=str)
    c_gene, t_gene = find(candidates, ["GENE", "GENES", "gene"]), find(truth, ["GENE", "gene"])
    if not c_gene or not t_gene: raise ValueError("Candidates and truth must contain GENE")
    c_id, t_id = find(candidates, ["SV_ID", "ID"]), find(truth, ["SV_ID", "ID"])
    candidates = candidates.reset_index(drop=True); candidates["_rank"] = candidates.index + 1
    detail = []
    for _, row in truth.iterrows():
        hits = candidates[candidates[c_gene].fillna("").str.upper() == str(row[t_gene]).upper()]
        if t_id and c_id and str(row.get(t_id, ".")) not in {"", ".", "nan"}:
            hits = hits[hits[c_id].astype(str) == str(row[t_id])]
        detail.append({"GENE": row[t_gene], "SV_ID": row.get(t_id, ".") if t_id else ".", "FOUND": "YES" if len(hits) else "NO", "BEST_RANK": int(hits["_rank"].min()) if len(hits) else "."})
    out = Path(a.output_tsv); out.parent.mkdir(parents=True, exist_ok=True); pd.DataFrame(detail).to_csv(out, sep="\t", index=False)
    found = sum(x["FOUND"] == "YES" for x in detail); summary = {"truth_count": len(detail), "found_count": found, "recall": found / len(detail) if detail else None}
    Path(a.output_json).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
