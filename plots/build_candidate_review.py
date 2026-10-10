#!/usr/bin/env python3
"""Create a small factual shortlist and a separate editable event-review table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from candidate_review import enrich_events, event_evidence, select_events, text
from plot_utils import read_tsv

REVIEW_COLUMNS = ["SAMPLE", "SV_ID", "GENE", "PANEL_STATUS", "STRUCTURAL_CHANGE", "GENE_RELEVANCE_CONTEXT", "REPORTED_CONTEXT", "SUPPORTING_OBSERVATIONS", "REVIEW_FLAGS", "WHAT_REMAINS_UNKNOWN", "WHY_GENE_MATTERS_REVIEW", "CONTRADICTORY_EVIDENCE_REVIEW", "NEXT_CHECK", "REVIEW_DECISION", "REVIEWER", "REVIEW_DATE", "SOURCE_CANDIDATE_TABLE", "SOURCE_EVIDENCE_TABLE"]


def write_review(events, sample, out, candidate_path, evidence_path, platform="lrs"):
    generated = []
    for _, row in events.iterrows():
        observations = event_evidence(row, platform=platform)
        context = [f"{d}: {label.replace(chr(10), ' ')}" for d, (state, label) in observations.items() if state == "REPORTED"]
        support = [f"{d}: {label.replace(chr(10), ' ')}" for d, (state, label) in observations.items() if state == "SUPPORTING"]
        flags = [f"{d}: {label.replace(chr(10), ' ')}" for d, (state, label) in observations.items() if state == "REVIEW"]
        unknown = [f"{d}: {label.replace(chr(10), ' ')}" for d, (state, label) in observations.items() if state in {"UNAVAILABLE", "UNRESOLVED"} or (d == "Population" and state == "NO_MATCH")]
        join_status = row.get("EVIDENCE_JOIN_STATUS", "NOT_REQUESTED")
        if join_status == "AMBIGUOUS_SV_GENE_ROWS":
            flags.append("Evidence join: ambiguous exact SV-gene source rows; not enriched")
        elif join_status == "NO_EXACT_SV_GENE_ROW":
            unknown.append("Evidence join: no exact SV-gene source row")
        record = {name: "" for name in REVIEW_COLUMNS}
        record.update({
            "SAMPLE": sample, "SV_ID": row["SV_ID"], "GENE": row["GENE"], "PANEL_STATUS": row["_panel_group"],
            "STRUCTURAL_CHANGE": f"{text(row, ['SVTYPE'])} {text(row, ['CHROM'])}:{text(row, ['START', 'POS'])}-{text(row, ['END'])}; {text(row, ['SV_GENE_EFFECT', 'SV_GENE_RELATIONSHIP'])}",
            "GENE_RELEVANCE_CONTEXT": text(row, ["FINAL_GENE_RELEVANCE_TIER", "GENE_RELEVANCE_TIER"]) + "; " + text(row, ["FINAL_PHENOTYPE_RANKING_SCOPE", "GENE_RELEVANCE_PHENOTYPE_SCOPE", "PHENOTYPE_SCORE_SCOPE"]),
            "REPORTED_CONTEXT": "; ".join(context), "SUPPORTING_OBSERVATIONS": "; ".join(support), "REVIEW_FLAGS": "; ".join(flags), "WHAT_REMAINS_UNKNOWN": "; ".join(unknown),
            "SOURCE_CANDIDATE_TABLE": str(Path(candidate_path).resolve()), "SOURCE_EVIDENCE_TABLE": str(Path(evidence_path).resolve()) if evidence_path else "",
        })
        generated.append(record)
    frame = pd.DataFrame(generated, columns=REVIEW_COLUMNS)
    frame.to_csv(out / f"{sample}_candidate_review.generated.tsv", sep="\t", index=False)
    editable = out / f"{sample}_candidate_review.tsv"
    # Never overwrite researcher notes. The generated table is refreshed
    # separately, so changes to a shortlist remain visible and reviewable.
    if not editable.exists():
        frame.to_csv(editable, sep="\t", index=False)
    else:
        print(f"[KEEP] existing editable review: {editable}")
    return frame


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--evidence-table")
    p.add_argument("--sample", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=5, help="Initial review rows per panel group")
    p.add_argument("--platform", choices=["lrs", "srs"], default="lrs")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    events = enrich_events(select_events(read_tsv(args.input), args.top_n), args.evidence_table)
    events.drop(columns=["_source_order"], errors="ignore").to_csv(out / f"{args.sample}_candidate_shortlist.tsv", sep="\t", index=False)
    records = write_review(events, args.sample, out, args.input, args.evidence_table, args.platform)
    read_evidence = "../../diagnostic_review/samplot" if args.platform == "srs" else "../13_breakpoint_evidence"
    context_limit = "prove that an overlapping repeat locus is expanded" if args.platform == "srs" else "prove that methylation is abnormal"
    guide = f"""START HERE: {args.sample}

1. Check coverage QC in ../01_qc and your existing alignment QC
   before reviewing candidates.
2. Open {args.sample}_candidate_review.tsv in Excel or a TSV viewer.
   It starts with up to {args.top_n} panel and {args.top_n} non-panel associations.
   Rows with unknown panel membership are shown separately if present.
   This is an initial workload choice, not a biological exclusion threshold.
3. Use {args.sample}_candidate_shortlist.tsv for copied pipeline fields.
   Track each review by SAMPLE + SV_ID + GENE, never by gene name alone.
4. Write WHY_GENE_MATTERS_REVIEW, CONTRADICTORY_EVIDENCE_REVIEW and NEXT_CHECK.
   Generated observations report computational/contextual findings. They do
   not establish pathogenicity or {context_limit}.
5. Compare the panel and non-panel matrices in ../05_candidate_prioritization.
   Availability is in evidence_availability/; mechanisms are in mechanisms/.
6. Inspect retained events in ../12_candidate_loci and the Samplot/read-level
   evidence in {read_evidence}. A gene BED alone has no exon detail.
7. Add lower-ranked direct disruptions or independent repeat/MEI findings
   when justified. This shortlist does not replace the complete results.

RERUNS: your editable review TSV is kept unchanged. The .generated.tsv file
and shortlist are refreshed. Reconcile new candidates by SAMPLE + SV_ID +
GENE; do not paste notes by row number. No selected candidate means an empty
shortlist, not proof that there is no causal event.

Candidate input: {Path(args.input).resolve()}
Detailed evidence: {Path(args.evidence_table).resolve() if args.evidence_table else 'Not supplied'}
Generated rows: {len(records)}
"""
    (out / "START_HERE.txt").write_text(guide, encoding="utf-8")
    print(f"[OK] review rows={len(records)}; start at {out / 'START_HERE.txt'}")


if __name__ == "__main__":
    main()
