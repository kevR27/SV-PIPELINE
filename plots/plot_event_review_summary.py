#!/usr/bin/env python3
"""Evidence availability and predicted mechanisms for the full candidate table."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from candidate_review import DOMAINS, domains_for, GROUP_LABELS, PANEL_GROUPS, STATE_COLORS, STATE_LABELS, evidence_table, normalise_events, text
from plot_utils import read_tsv, save_figure, set_thesis_style
from candidate_review import event_frequency_table, MISSING
from plot_caller_concordance import plot_upset

EVENT_DOMAINS = {"Call support", "Population", "Local depth", "SV phase", "Straglr", "TLDR", "Methylation", "CNVpytor depth", "GRIDSS", "Repeat locus", "MELT"}


def evidence_membership(raw, shortlist, gene_bed=None):
    """Describe overlapping evidence sets in the complete master table.

    Any gene association can supply a panel flag. HPO evidence counts only when
    its scope explicitly says patient-specific; cohort seed terms do not qualify.
    """
    selected = set(shortlist["SV_ID"].astype(str))
    frequencies = event_frequency_table(raw).set_index("SV_ID")
    gene_index = None
    if gene_bed:
        from plot_genome_wide_sv_map import chrom_name, interval_index
        bed = pd.read_csv(gene_bed,sep="\t",comment="#",header=None,usecols=[0,1,2],names=["CHROM","START0","END0"])
        bed["CHROM"] = bed["CHROM"].map(chrom_name)
        gene_index = interval_index(bed)
    rows = []
    for sv_id, group in raw.groupby("SV_ID", sort=False):
        if str(sv_id).upper() in MISSING:
            raise ValueError("Counts need nonmissing master SV IDs")
        caller = pd.to_numeric(group.get("CALLER_COUNT", pd.Series(dtype=float)), errors="coerce").dropna().unique()
        multi = bool(caller[0] >= 2) if len(caller) == 1 else pd.NA
        genic = []
        panel = []
        phenotype = []
        for row in group.to_dict("records"):
            gene = text(row, ["GENE", "GENES", "ANNotsv_Gene", "Gene"])
            effect = text(row, ["SV_GENE_EFFECT", "SV_GENE_RELATIONSHIP"]).upper()
            # An annotated nearby gene is not necessarily an overlapping gene.
            explicit_context = effect != "." and not any(x in effect for x in ["PROXIMAL", "NEAR_GENE", "INTERGENIC", "UNRESOLVED", "NO_GENE"])
            genic.append(pd.NA if effect == "." and gene != "." else gene != "." and explicit_context)
            membership = text(row, ["PANEL_STATUS", "panel_gene"]).upper()
            panel.append(True if membership in {"PANEL_GENE","YES","TRUE","1"} else False if membership in {"NONPANEL_GENE","NON_PANEL","NO","FALSE","0"} else pd.NA)
            scope = text(row, ["GENE_RELEVANCE_PHENOTYPE_SCOPE", "FINAL_PHENOTYPE_RANKING_SCOPE", "PHENOTYPE_SCORE_SCOPE"]).upper()
            score = pd.to_numeric(text(row, ["PATIENT_HPO_OVERLAP", "PATIENT_HPO_SCORE", "PHENOTYPE_SCORE"]), errors="coerce")
            phenotype.append(bool(score > 0) if "PATIENT" in scope and "NO_PATIENT" not in scope and "UNAVAILABLE" not in scope else pd.NA)
        genic_flag = True if any(x is True for x in genic) else pd.NA if any(x is pd.NA for x in genic) else False
        if gene_index is not None:
            from plot_genome_wide_sv_map import chrom_name, overlaps
            geometry = [col for col in ["CHROM","START","END","SVTYPE","CHR2","POS2"] if col in group]
            row = group.iloc[0]
            start = pd.to_numeric(row.get("START"),errors="coerce")
            end = pd.to_numeric(row.get("END"),errors="coerce")
            chrom = chrom_name(row.get("CHROM"))
            genic_flag = pd.NA
            if pd.notna(start) and start >= 1 and not group[geometry].nunique(dropna=False).gt(1).any():
                kind = str(row.get("SVTYPE")).upper()
                right = end if kind in {"DEL","DUP","INV","CNV"} and pd.notna(end) and end >= start else start
                genic_flag = bool(overlaps(gene_index,chrom,int(start)-1,int(right)))
                if kind in {"BND","TRA","CTX"}:
                    pos2 = pd.to_numeric(row.get("POS2"),errors="coerce")
                    if pd.notna(pos2) and pos2 >= 1:
                        genic_flag |= bool(overlaps(gene_index,chrom_name(row.get("CHR2")),int(pos2)-1,int(pos2)))
        patient_flag = True if any(x is True for x in phenotype) else pd.NA if any(x is pd.NA for x in phenotype) else False
        panel_flag = True if any(x is True for x in panel) else pd.NA if any(x is pd.NA for x in panel) else False
        panel_or_hpo = True if panel_flag is True or patient_flag is True else pd.NA if panel_flag is pd.NA or patient_flag is pd.NA else False
        af_status = frequencies.at[sv_id, "AF_STATUS"]
        rows.append({"SV_ID": sv_id, "Multi-caller": multi,
                     "AF ≤1% or unknown": af_status in {"AF_UNKNOWN", "MEASURED_AT_OR_BELOW_CUTOFF"},
                     "Gene overlap": genic_flag, "Panel": panel_flag, "Patient HPO": patient_flag, "Panel or patient HPO": panel_or_hpo,
                     "Review shortlist": str(sv_id) in selected, "AF_STATUS": af_status})
    return pd.DataFrame(rows, columns=["SV_ID", "Multi-caller", "AF ≤1% or unknown", "Gene overlap", "Panel", "Patient HPO", "Panel or patient HPO", "Review shortlist", "AF_STATUS"])


def plot_evidence_counts(raw, shortlist, out, sample, caller_summary=None, raw_vcfs=(), gene_bed=None):
    out.mkdir(parents=True, exist_ok=True)
    membership = evidence_membership(raw, shortlist, gene_bed)
    membership.to_csv(out / f"{sample}_evidence_membership.tsv", sep="\t", index=False)
    counts = []
    raw_total = 0
    for item in raw_vcfs:
        name, path = item.split("=", 1)
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as handle:
            n = sum(bool(line.strip()) and not line.startswith("#") for line in handle)
        raw_total += n
        counts.append({"GROUP": "Raw "+name, "COUNT": n, "UNEVALUABLE": 0,
                       "UNIT": "CALLER_VCF_RECORDS", "SCOPE": "Includes duplicate reports and separate breakend records", "SOURCE": path})
    counts.insert(0, {"GROUP": "Raw caller records (sum)", "COUNT": raw_total if raw_vcfs else np.nan,
                      "UNEVALUABLE": 0, "UNIT": "CALLER_VCF_RECORDS", "SCOPE": "Not comparable to unique merged SVs", "SOURCE": ";".join(raw_vcfs) or "Unavailable"})
    if caller_summary:
        merged = read_tsv(caller_summary)
        counts.append({"GROUP": "Merged breakpoint SVs", "COUNT": merged["SV_ID"].nunique(), "UNEVALUABLE": 0,
                       "UNIT": "UNIQUE_MASTER_SV", "SCOPE": "Breakpoint caller union before depth-only additions", "SOURCE": str(caller_summary)})
    counts.append({"GROUP": "Annotated master SVs", "COUNT": len(membership), "UNEVALUABLE": 0,
                   "UNIT": "UNIQUE_MASTER_SV", "SCOPE": "Complete supplied table; SRS may include depth-only events", "SOURCE": "Complete evidence table"})
    for column in ["Multi-caller", "AF ≤1% or unknown", "Gene overlap", "Panel", "Patient HPO", "Panel or patient HPO", "Review shortlist"]:
        counts.append({"GROUP": column, "COUNT": int(membership[column].eq(True).sum()),
                       "UNEVALUABLE": int(membership[column].isna().sum()), "UNIT": "UNIQUE_MASTER_SV",
                       "SCOPE": "Independent overlapping evidence group, not a sequential filter", "SOURCE": "Complete evidence table / existing shortlist"})
    frame = pd.DataFrame(counts)
    frame.to_csv(out / f"{sample}_counts_at_each_step.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1, 2]})
    for ax, unit, title in zip(axes, ["CALLER_VCF_RECORDS", "UNIQUE_MASTER_SV"], ["Before merging: caller records", "After merging: unique master SVs"]):
        sub = frame[frame["UNIT"].eq(unit)].reset_index(drop=True)
        ax.barh(sub["GROUP"], sub["COUNT"].fillna(0), color="#999999" if unit == "CALLER_VCF_RECORDS" else "#0072B2", height=.65)
        for i, row in sub.iterrows():
            label = "Unavailable" if pd.isna(row["COUNT"]) else f"{int(row['COUNT']):,}"
            if row["UNEVALUABLE"]:
                label += f"; {int(row['UNEVALUABLE']):,} unknown"
            ax.text(0 if pd.isna(row["COUNT"]) else row["COUNT"], i, "  "+label, va="center", fontsize=9)
        ax.invert_yaxis()
        ax.set_xlim(0, max(sub["COUNT"].fillna(0).max(), 1)*1.5)
        ax.set_title(title, loc="left", fontsize=12)
        ax.set_xlabel("Records" if unit == "CALLER_VCF_RECORDS" else "Unique SVs")
    fig.suptitle("Counts at each analysis stage / evidence group", fontweight="bold")
    fig.text(.5,.02,"The right-hand bars are overlapping groups, not a sequential attrition funnel. The shortlist is an adjustable review workload, not validated diagnoses.\nAF ≤1% and unknown AF are retained together. Patient HPO needs explicit patient-specific evidence; cohort seed terms do not count.",ha="center",fontsize=9)
    fig.tight_layout(rect=[0,.10,1,.94])
    save_figure(fig, out / f"{sample}_counts_at_each_step")
    plt.close(fig)
    # A nonsequential alternative: exact intersections of evaluated evidence.
    layers = ["Multi-caller", "AF ≤1% or unknown", "Gene overlap", "Panel", "Review shortlist"]
    evaluable = membership.dropna(subset=layers)
    vectors = evaluable[layers].astype(int).astype(str).agg("".join, axis=1).value_counts()
    intersections = pd.DataFrame([{"SUPP_VEC": vec, "intersection_size": n, "caller_count": vec.count("1")} for vec,n in vectors.items()])
    if not intersections.empty:
        intersections = intersections.sort_values("intersection_size",ascending=False)
        intersections.to_csv(out / f"{sample}_evidence_intersections.tsv",sep="\t",index=False)
        sizes = pd.DataFrame({"caller": layers, "size": [int(evaluable[col].sum()) for col in layers]})
        fig = plot_upset(intersections, sizes, layers, f"Overlapping SV evidence groups ({len(evaluable)} fully evaluable SVs)")
        fig.axes[1].set_ylabel("Unique master SVs")
        fig.axes[3].set_yticklabels([])  # The same group names are already shown beside the set sizes.
        fig.axes[3].tick_params(axis="y",length=0)
        fig.axes[3].set_xticklabels(range(1,len(intersections)+1),rotation=0)
        fig.texts[-1].set_text(f"Exact intersections of independently defined evidence groups. {len(membership)-len(evaluable)} SVs omitted because at least one group is unknown.\nColumn numbers match the rows of the source intersection TSV; no sequential filtering is implied.")
        save_figure(fig, out / f"{sample}_evidence_upset")
        plt.close(fig)


def availability_counts(observations, domains=None):
    rows = []
    for domain in DOMAINS if domains is None else domains:
        sub = observations[observations["DOMAIN"].eq(domain)].copy()
        unit = "UNIQUE_SV" if domain in EVENT_DOMAINS else "SV_GENE_ASSOCIATION"
        if unit == "UNIQUE_SV":
            # Event-wide evidence must agree across gene rows. Conflicting
            # propagated observations stay unresolved rather than picking one.
            sub = sub.groupby("SV_ID", as_index=False).agg(STATE=("STATE", lambda x: x.iloc[0] if x.nunique() == 1 else "UNRESOLVED"))
        total = len(sub)
        for state in STATE_COLORS:
            count = int(sub["STATE"].eq(state).sum())
            rows.append({"DOMAIN": domain, "COUNTING_UNIT": unit, "STATE": state, "COUNT": count, "TOTAL": total, "PERCENT": 100 * count / total if total else 0})
    return pd.DataFrame(rows)


def mechanism_category(row):
    effect = text(row, ["SV_GENE_EFFECT", "SV_GENE_RELATIONSHIP"]).upper()
    svtype = text(row, ["SVTYPE"]).upper()
    if effect in {"WHOLE_GENE_DELETION", "WHOLE_GENE_DOSAGE_CONTEXT"} and svtype == "DEL":
        return "Whole-gene loss"
    if effect in {"WHOLE_GENE_DUPLICATION", "WHOLE_GENE_DOSAGE_CONTEXT"} and svtype == "DUP":
        return "Whole-gene copy gain"
    if "PARTIAL_GENE" in effect:
        return "Partial gene loss" if svtype == "DEL" else "Partial gene duplication" if svtype == "DUP" else "Partial gene overlap"
    if "SPANNED_BY_INVERSION" in effect or "SPANS_INTACT_GENE" in effect:
        return "Gene spanned by inversion"
    if "BREAKPOINT_WITHIN" in effect or "BREAKPOINT_IN_GENE" in effect:
        return "Breakpoint within gene"
    if "INSERTION_WITHIN" in effect or effect == "INSERTION_IN_GENE":
        return "Insertion within gene"
    if "INTRON" in effect or text(row, ["VEP_TRANSCRIPT_REGION_CLASS"]).upper() == "INTRONIC_ONLY":
        return "Intronic context"
    if "PROXIMAL" in effect or "NEAR_GENE" in effect:
        return "Proximal context"
    if "CONTEXT" in effect or effect == "GENE_FULLY_SPANNED":
        return "Other interval context"
    return "Unresolved / other"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--platform", choices=["lrs", "srs"], default="lrs")
    p.add_argument("--master-table", help="Complete integrated table for independent counts")
    p.add_argument("--selection-table", help="Existing review shortlist; defines the plotted review group")
    p.add_argument("--caller-summary")
    p.add_argument("--gene-bed", help="Actual gene intervals, for coordinate-based overlap counts")
    p.add_argument("--raw-vcf", action="append", default=[], metavar="CALLER=PATH")
    args = p.parse_args()
    set_thesis_style()
    work = normalise_events(read_tsv(args.input))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if args.master_table and args.selection_table:
        plot_evidence_counts(read_tsv(args.master_table), read_tsv(args.selection_table), out / "analysis_counts", args.sample, args.caller_summary, args.raw_vcf, args.gene_bed)
    availability_out = out / "evidence_availability"
    mechanism_out = out / "mechanisms"
    availability_out.mkdir(exist_ok=True)
    mechanism_out.mkdir(exist_ok=True)
    domains = domains_for(args.platform)
    counts = availability_counts(evidence_table(work, platform=args.platform), domains)
    counts.to_csv(availability_out / f"{args.sample}_evidence_availability.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(13, 8.5))
    left = np.zeros(len(domains))
    for state in STATE_COLORS:
        values = counts[counts["STATE"].eq(state)].set_index("DOMAIN")["PERCENT"].reindex(domains).to_numpy()
        ax.barh(domains, values, left=left, label=STATE_LABELS[state], color=STATE_COLORS[state], edgecolor="white", height=0.7)
        left += values
    labels = [d + (" [SV]" if d in EVENT_DOMAINS else " [SV–gene]") for d in domains]
    ax.set_yticks(range(len(domains)), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percentage of the stated counting unit")
    ax.set_title(f"{args.sample}: evidence availability in the candidate table", loc="left", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.22, -0.15), ncol=3, frameon=False, fontsize=9)
    fig.text(0.5, 0.012, "SV-wide domains count each SV once. Gene-specific domains count SV–gene associations.\nUnavailable means absent or unevaluable in this table; it is not evidence against the event. Conflicting event-wide rows are unresolved.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.23, right=0.97, top=0.92, bottom=0.26)
    save_figure(fig, availability_out / f"{args.sample}_evidence_availability", dpi=300)
    plt.close(fig)

    work["PREDICTED_MECHANISM"] = [mechanism_category(row) for _, row in work.iterrows()]
    summary = work.groupby(["_panel_group", "PREDICTED_MECHANISM"], as_index=False).agg(SV_GENE_ASSOCIATIONS=("SV_ID", "size"), UNIQUE_SVS=("SV_ID", "nunique"), UNIQUE_GENES=("GENE", "nunique"))
    summary = summary.rename(columns={"_panel_group": "PANEL_STATUS"})
    summary.to_csv(mechanism_out / f"{args.sample}_mechanism_summary.tsv", sep="\t", index=False)
    categories = list(dict.fromkeys(summary.sort_values("SV_GENE_ASSOCIATIONS", ascending=False)["PREDICTED_MECHANISM"]))
    fig, ax = plt.subplots(figsize=(13, max(5.3, len(categories)*0.48+2)))
    left = np.zeros(len(categories))
    for group, color in zip(PANEL_GROUPS, ["#0072B2", "#E69F00", "#999999"]):
        sub = summary[summary["PANEL_STATUS"].eq(group)].set_index("PREDICTED_MECHANISM")
        values = sub["SV_GENE_ASSOCIATIONS"].reindex(categories, fill_value=0).to_numpy()
        ax.barh(categories, values, left=left, label=GROUP_LABELS[group], color=color, height=0.7)
        left += values
    for i, value in enumerate(left):
        ax.text(value, i, f"  {int(value):,}", va="center", fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, max(left, default=0)*1.15+1)
    ax.set_xlabel("Unique SV–gene associations")
    ax.set_title(f"{args.sample}: predicted SV–gene mechanisms", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right")
    fig.text(0.5, 0.012, "Categories describe annotated geometry, not demonstrated loss of function. One SV may affect several genes and contribute to several categories.\nPartial gene loss does not establish coding-exon loss. Unique SV counts per category are provided separately in the TSV and are not additive.", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.29, right=0.97, top=0.9, bottom=0.19)
    save_figure(fig, mechanism_out / f"{args.sample}_mechanism_summary", dpi=300)
    plt.close(fig)
    print(f"[OK] {len(work)} unique SV-gene associations; review summaries: {out}")


if __name__ == "__main__":
    main()
