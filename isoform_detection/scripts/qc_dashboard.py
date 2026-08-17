#!/usr/bin/env python3
"""Builds the QC dashboard: one HTML page gathering the per-cell-type QC
reports into a single picture.
 
The page has, from top to bottom:
  - a summary panel of the headline cohort numbers
  - tag completeness by cell type (bar)
  - completeness vs reads per cell type (scatter) with the correlation
  - parameter sensitivity: mean retention vs ZW cutoff (line)
  - read fate by cell type (stacked bar) -- only shown when IC:i:0 profiling
    data is present, since that is an optional pipeline step
 
Each chart carries a short caption stating what it counts, because the charts use
different denominators and would otherwise appear to contradict one another.
"""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import argparse
import glob
import os
import math
from datetime import datetime


# Small helper so each chart can carry a plain-English caption explaining what it is actually counting.
def caption(text):
    return (
        "<p style='max-width: 900px; margin: 0 auto 30px auto; "
        "color: #555; font-size: 14px; text-align: center;'>"
        + text +
        "</p>"
    )

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return None
    return cov / (sx * sy)

def main():
    parser = argparse.ArgumentParser(description="Generate the master QC dashboard")
    parser.add_argument("-c", "--comp_dir", required=True, help="Directory with tag-completeness TSVs")
    parser.add_argument("-s", "--sens_dir", required=True, help="Directory with sensitivity TSVs")
    parser.add_argument("-o", "--output", required=True, help="Output HTML file path")
    parser.add_argument("-i", "--ic0_dir", default=None,
                        help="Optional directory with IC:i:0 profiling TSVs (for the read-fate chart)")
    parser.add_argument("-t", "--threshold", type=float, default=0.8,
                        help="ZW assignment threshold, used only for chart labels [%(default)s]")
    parser.add_argument("-r", "--run-identity", default=None,
                        help="Run name, shown on the page so an exported chart can be traced back")
    args = parser.parse_args()
 
    html_parts = []

# Load completeness
    comp_files = glob.glob(os.path.join(args.comp_dir, "*_completeness.tsv"))
    df_comp = None
    if comp_files:
        df_comp = pd.concat([pd.read_csv(f, sep="\t") for f in comp_files], ignore_index=True)
        df_comp = df_comp.sort_values(by="Sample")
 
# Load sensitivity
    sens_files = glob.glob(os.path.join(args.sens_dir, "*_sensitivity.tsv"))
    df_sens = None
    if sens_files:
        df_sens = pd.concat([pd.read_csv(f, sep="\t") for f in sens_files], ignore_index=True)
 

# SUMMARY PANEL (headline numbers)
    if df_comp is not None:
        total = int(df_comp["Total_Reads"].sum())
        assigned = int(df_comp["Assigned_Reads"].sum())
        pooled = assigned / total * 100 if total else 0
        log_reads = [math.log10(d) for d in df_comp["Total_Reads"] if d > 0]
        compl = [c for d, c in zip(df_comp["Total_Reads"], df_comp["Completeness_Percentage"]) if d > 0]
        r_reads = pearson(log_reads, compl)
 
        swing_txt = ""
        if df_sens is not None:
            means = df_sens.groupby("Threshold")["Retention_Rate"].mean()
            lo_t, hi_t = means.index[0], means.index[-1]
            swing = means.iloc[0] - means.iloc[-1]
            swing_txt = f"<b>{swing:.1f} pts</b> retention swing ({lo_t:g}&ndash;{hi_t:g})"
 
        def stat(label, value):
            return (
                "<div style='display:inline-block; text-align:center; margin:0 25px;'>"
                f"<div style='font-size:26px; font-weight:bold; color:#1f3864;'>{value}</div>"
                f"<div style='font-size:13px; color:#666;'>{label}</div>"
                "</div>"
            )
 
        panel = "<div style='text-align:center; margin:20px auto 35px auto; padding:18px; " \
                "background:#f5f7fb; border-radius:8px; max-width:900px;'>"
        panel += stat("cell types", len(df_comp))
        panel += stat("reads assigned", f"{pooled:.1f}%")
        panel += stat("total reads", f"{total/1e6:.0f}M")
        if r_reads is not None:
            panel += stat("reads correlation", f"{r_reads:.2f}")
        panel += "</div>"
        if swing_txt:
            panel += f"<p style='text-align:center; color:#555; margin-top:-20px;'>{swing_txt}</p>"
        html_parts.append(panel)
 

# CHART 1: tag completeness by cell type
    if df_comp is not None:
        fig = px.bar(
            df_comp, x="Sample", y=["Assigned_Reads", "Unassigned_Reads"],
            title="Tag completeness by cell type",
            color_discrete_map={"Assigned_Reads": "#2ca02c", "Unassigned_Reads": "#d62728"},
        )
        html_parts.append(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        html_parts.append(caption(
            "Each bar shows all reads in one cell type, split into assigned (green) and "
            "unassigned (red). Denominator: all reads in the cell type."
        ))
        html_parts.append("<hr>")
 

# CHART 2: completeness by cell type
    if df_comp is not None:
        d = df_comp.copy()
        r_txt = f" (r = {r_reads:.2f} on log reads)" if r_reads is not None else ""
        fig = px.scatter(
            d, x="Total_Reads", y="Completeness_Percentage",
            color="Sample", text="Sample", log_x=True,
            title="Completeness vs total reads per cell type" + r_txt,
        )
        fig.update_traces(textposition="top center", textfont_size=9, marker_size=11)
        fig.update_layout(xaxis_title="Total reads per cell type (log scale)",
                          yaxis_title="Tag completeness (%)", showlegend=False)
        html_parts.append(fig.to_html(full_html=False, include_plotlyjs=False))
        html_parts.append(caption(
            "Each point is one cell type. The trend shows how assignment rate rises with "
            "the total number of reads assigned to that cell type. Note this is a count of pooled reads, not a per-cell sequencing depth: a cell type with many cells at low coverage and one with few cells at high coverage can carry the same total."
        ))
        html_parts.append("<hr>")
 

# CHART 3: parameter sensitivity
    if df_sens is not None:
        df_sens_avg = df_sens.groupby("Threshold")["Retention_Rate"].mean().reset_index()
        fig = px.line(
            df_sens_avg, x="Threshold", y="Retention_Rate", markers=True,
            title="Parameter sensitivity: average read retention vs ZW cutoff",
            labels={"Retention_Rate": "Avg reads retained (%)",
                    "Threshold": "kallisto ZW probability cutoff"},
        )
        fig.update_traces(line_color="#1f77b4", marker=dict(size=10))
        html_parts.append(fig.to_html(full_html=False, include_plotlyjs=False))
        html_parts.append(caption(
            "How many reads survive as the confidence cutoff is raised, averaged across all "
            "cell types. A flat line means robustness to the cutoff. Denominator: "
            "kallisto-aligned reads only, which is why this differs from the completeness chart."
        ))
        html_parts.append("<hr>")
 

# CHART 4: read fate by cell type (only if IC:i:0 data present)
    ic0_df = None
    if args.ic0_dir and os.path.isdir(args.ic0_dir):
        ic0_files = glob.glob(os.path.join(args.ic0_dir, "*_ic0.tsv"))
        all_file = os.path.join(args.ic0_dir, "ALL_ic0.tsv")
        if os.path.exists(all_file):
            ic0_df = pd.read_csv(all_file, sep="\t")
        elif ic0_files:
            ic0_df = pd.concat([pd.read_csv(f, sep="\t") for f in ic0_files], ignore_index=True)
 
    if ic0_df is not None and not ic0_df.empty:
        d = ic0_df.copy()
        d["tot"] = d["Assigned"] + d["Scored_But_Low"] + d["Not_Pseudoaligned"]
        d["p_assigned"] = 100 * d["Assigned"] / d["tot"]
        d["p_low"] = 100 * d["Scored_But_Low"] / d["tot"]
        d["p_np"] = 100 * d["Not_Pseudoaligned"] / d["tot"]
        d = d.sort_values("tot", ascending=False)
 
        fig = go.Figure()
        fig.add_bar(x=d["Sample"], y=d["p_assigned"], name="Assigned (IC > 0)", marker_color="#2ca02c")
        fig.add_bar(x=d["Sample"], y=d["p_low"], name=f"Unassigned: scored but below {args.threshold:g}", marker_color="#ff7f0e")
        fig.add_bar(x=d["Sample"], y=d["p_np"], name="Unassigned: never pseudo-aligned", marker_color="#d62728")
        fig.update_layout(
            barmode="stack",
            title="Read fate by cell type",
            yaxis_title="Percentage of reads (%)",
            xaxis_title=None,
            legend=dict(orientation="h", yanchor="bottom", y=-0.35),
        )
        html_parts.append(fig.to_html(full_html=False, include_plotlyjs=False))
        html_parts.append(caption(
            "Every read in each cell type, split three ways: assigned (green), unassigned but "
            f"scored below {args.threshold:g} (orange), and unassigned because kallisto never pseudo-aligned it "
            "(red). The red fraction shows reads that matched no assembled transcript."
        ))
    else:
        html_parts.append(caption(
            "Read-fate chart not shown: it requires IC:i:0 profiling data. Run the pipeline "
            "with run_ic0_profiling enabled, and pass the ic0 directory to this script, to see it."
        ))
 
  
  # WRITE PAGE
    with open(args.output, "w") as f:
        page_title = "Alternative Splicing QC Dashboard"
        if args.run_identity:
            page_title += f" \u2014 {args.run_identity}"
        f.write(f"<html><head><title>{page_title}</title></head>")
        f.write("<body style='font-family: Verdana, sans-serif;'>")
        f.write("<h1 style='text-align: center; padding-top: 20px; margin-bottom: 4px;'>"
                "Alternative Splicing Pipeline - Quality Control Dashboard</h1>")
        # Run name and timestamp, so a chart exported from this page can be traced
        # back to the run that produced it.
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
        subtitle = f"generated {stamp}"
        if args.run_identity:
            subtitle = f"run: {args.run_identity} &nbsp;|&nbsp; " + subtitle
        f.write("<p style='text-align: center; color: #666; font-size: 14px; "
                f"margin-top: 0; margin-bottom: 24px;'>{subtitle}</p>")
        for part in html_parts:
            f.write(part)
        f.write("</body></html>")
 
    print(f"QC dashboard successfully generated at: {args.output}")
 
 
if __name__ == "__main__":
    main()