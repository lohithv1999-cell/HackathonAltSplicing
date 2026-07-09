#!/usr/bin/env python3
"""QC summary statistics.

Reads all the per-cell-type completeness and sensitivity TSVs and computes the
cohort-level headline numbers: the pooled assignment rate, the correlation between
sequencing depth and completeness, and the average sensitivity retention at each
threshold. Writes them to one summary TSV so every run documents its own figures.
"""

import argparse
import csv
import glob
import math
import os


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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
    parser = argparse.ArgumentParser(description="Compute cohort-level QC summary statistics")
    parser.add_argument("-c", "--comp_dir", required=True, help="Directory with tag-completeness TSVs")
    parser.add_argument("-s", "--sens_dir", required=True, help="Directory with sensitivity TSVs")
    parser.add_argument("-o", "--output", required=True, help="Output summary TSV")
    args = parser.parse_args()

    comp_rows = []
    for f in glob.glob(os.path.join(args.comp_dir, "*_completeness.tsv")):
        comp_rows.extend(read_tsv(f))

    summary = []

    if comp_rows:
        total = sum(int(r["Total_Reads"]) for r in comp_rows)
        assigned = sum(int(r["Assigned_Reads"]) for r in comp_rows)
        pooled_rate = (assigned / total * 100) if total > 0 else 0.0

        depths = [int(r["Total_Reads"]) for r in comp_rows if int(r["Total_Reads"]) > 0]
        compl = [float(r["Completeness_Percentage"]) for r in comp_rows if int(r["Total_Reads"]) > 0]
        r_depth = pearson([math.log10(d) for d in depths], compl)

        pcts = [float(r["Completeness_Percentage"]) for r in comp_rows]
        names = [r["Sample"] for r in comp_rows]
        i_min = min(range(len(pcts)), key=lambda i: pcts[i])
        i_max = max(range(len(pcts)), key=lambda i: pcts[i])

        summary += [
            ("n_cell_types", len(comp_rows)),
            ("total_reads", total),
            ("total_assigned", assigned),
            ("pooled_assignment_rate_pct", round(pooled_rate, 2)),
            ("completeness_corr_log_depth", round(r_depth, 3) if r_depth is not None else "NA"),
            ("completeness_min_pct", f"{pcts[i_min]:.2f} ({names[i_min]})"),
            ("completeness_max_pct", f"{pcts[i_max]:.2f} ({names[i_max]})"),
        ]

    sens_rows = []
    for f in glob.glob(os.path.join(args.sens_dir, "*_sensitivity.tsv")):
        sens_rows.extend(read_tsv(f))

    if sens_rows:
        by_threshold = {}
        for r in sens_rows:
            t = float(r["Threshold"])
            by_threshold.setdefault(t, []).append(float(r["Retention_Rate"]))
        for t in sorted(by_threshold):
            vals = by_threshold[t]
            mean_ret = sum(vals) / len(vals)
            summary.append((f"mean_retention_at_{t:g}_pct", round(mean_ret, 2)))

        ts = sorted(by_threshold)
        if len(ts) >= 2:
            lo_mean = sum(by_threshold[ts[0]]) / len(by_threshold[ts[0]])
            hi_mean = sum(by_threshold[ts[-1]]) / len(by_threshold[ts[-1]])
            summary.append(("retention_swing_pct", round(lo_mean - hi_mean, 2)))

    with open(args.output, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Metric", "Value"])
        for metric, value in summary:
            w.writerow([metric, value])

    print("QC summary statistics:")
    for metric, value in summary:
        print(f"  {metric}: {value}")


if __name__ == "__main__":
    main()