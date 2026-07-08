#!/usr/bin/env python3
"""Builds the QC dashboard: one HTML page with two interactive charts that
gather the per-cell-type QC reports into a single picture.

Chart 1 (Tag completeness): for each cell type, how many reads got a real isoform
tag versus how many did not. The percentage is out of ALL reads in that cell type.

Chart 2 (Parameter sensitivity): how the average read-retention rate changes as the
ZW confidence cutoff is raised. This percentage is out of KALLISTO-ALIGNED reads
only, which is a different denominator from chart 1 (so the two charts are not
expected to show the same number at 0.8). A short caption under each chart spells
this out so the difference is clear to anyone reading the page.
"""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import argparse
import glob
import os


# Small helper so each chart can carry a plain-English caption explaining what it
# is actually counting.
def caption(text):
    return (
        "<p style='max-width: 900px; margin: 0 auto 30px auto; "
        "color: #555; font-size: 14px; text-align: center;'>"
        + text +
        "</p>"
    )


def main():
    parser = argparse.ArgumentParser(description="Generate the master QC dashboard")
    parser.add_argument("-c", "--comp_dir", required=True, help="Directory with tag-completeness TSVs")
    parser.add_argument("-s", "--sens_dir", required=True, help="Directory with sensitivity TSVs")
    parser.add_argument("-o", "--output", required=True, help="Output HTML file path")
    args = parser.parse_args()

    # Chart 1: tag completeness
    comp_files = glob.glob(os.path.join(args.comp_dir, "*_completeness.tsv"))
    if comp_files:
        df_comp = pd.concat([pd.read_csv(f, sep="\t") for f in comp_files], ignore_index=True)
        df_comp = df_comp.sort_values(by="Sample")
        fig_comp = px.bar(
            df_comp,
            x="Sample",
            y=["Assigned_Reads", "Unassigned_Reads"],
            title="Tag completeness by cell type (default cutoff 0.8)",
            color_discrete_map={"Assigned_Reads": "#2ca02c", "Unassigned_Reads": "#d62728"},
        )
        comp_caption = caption(
            "Each bar shows all reads in one cell type, split into those that received "
            "a real isoform tag (green) and those that did not (red). The percentage of "
            "green is the 'completeness' for that cell type. Denominator: all reads in "
            "the cell type."
        )
    else:
        fig_comp = go.Figure().update_layout(title="No completeness data found")
        comp_caption = ""

    # Chart 2: parameter sensitivity
    sens_files = glob.glob(os.path.join(args.sens_dir, "*_sensitivity.tsv"))
    if sens_files:
        df_sens = pd.concat([pd.read_csv(f, sep="\t") for f in sens_files], ignore_index=True)
        # Average the per-cell-type retention rate at each threshold.
        df_sens_avg = df_sens.groupby("Threshold")["Retention_Rate"].mean().reset_index()
        fig_sens = px.line(
            df_sens_avg,
            x="Threshold",
            y="Retention_Rate",
            markers=True,
            title="Parameter sensitivity: average read retention vs ZW cutoff",
            labels={
                "Retention_Rate": "Avg reads retained (%)",
                "Threshold": "kallisto ZW probability cutoff",
            },
        )
        fig_sens.update_traces(line_color="#1f77b4", marker=dict(size=10))
        sens_caption = caption(
            "This line shows how many reads survive as the confidence cutoff is raised, "
            "averaged across all cell types. A flat line means the pipeline is robust to "
            "the cutoff choice. Denominator: kallisto-aligned reads only, which is why "
            "this percentage differs from the completeness chart above."
        )
    else:
        fig_sens = go.Figure().update_layout(title="No sensitivity data found")
        sens_caption = ""

    # Write the page
    with open(args.output, "w") as f:
        f.write("<html><head><title>Alternative Splicing QC Dashboard</title></head>")
        f.write("<body style='font-family: Verdana, sans-serif;'>")
        f.write("<h1 style='text-align: center; padding-top: 20px;'>"
                "Alternative Splicing Pipeline - Quality Control Dashboard</h1>")

        f.write(fig_comp.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(comp_caption)

        f.write("<hr>")

        f.write(fig_sens.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(sens_caption)

        f.write("</body></html>")

    print(f"QC dashboard successfully generated at: {args.output}")


if __name__ == "__main__":
    main()