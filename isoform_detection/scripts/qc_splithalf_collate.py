#!/usr/bin/env python3
"""Gather the per-cell-type split-half stability results into one table.

Sorted by read count, since that is the obvious thing to check the figures against.
"""

import argparse
import csv
import os


def main():
    parser = argparse.ArgumentParser(description="Collate split-half stability TSVs")
    parser.add_argument("-i", "--inputs", nargs="+", required=True, help="Per-cell-type stability TSVs")
    parser.add_argument("-o", "--output", required=True, help="Combined output TSV")
    args = parser.parse_args()

    rows = []
    fieldnames = None
    for path in args.inputs:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if fieldnames is None:
                fieldnames = reader.fieldnames
            rows.extend(reader)

    if not rows:
        raise SystemExit("No stability rows found")

    def read_total(row):
        try:
            return float(row.get("Total_Counts_A", 0)) + float(row.get("Total_Counts_B", 0))
        except (TypeError, ValueError):
            return 0.0

    rows.sort(key=depth, reverse=True)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        recovered = row.get("Pct_MultiExon_A_Recovered", row.get("Pct_A_Recovered_In_B", "NA"))

    # Cell types the compare step could not assess at all, usually because the
    # assembly had no multi-exon transcripts. Kept out of the averages.
    skipped = [r for r in rows if r.get("Status", "ok") != "ok"]
    if skipped:
        print(f"\n  {len(skipped)} of {len(rows)} cell types could not be assessed:")
        reasons = {}
        for r in skipped:
            reasons.setdefault(r["Status"], []).append(r["Sample"])
        for reason, samples in reasons.items():
            shown = ", ".join(samples[:5])
            more = f" and {len(samples) - 5} others" if len(samples) > 5 else ""
            print(f"    {reason}: {shown}{more}")

    usable = [float(r["Counts_Pearson"]) for r in rows
              if r.get("Counts_Pearson") not in (None, "", "NA")]
    if usable:
        print(f"\n  mean Pearson across {len(usable)} cell types: {sum(usable) / len(usable):.4f}")
        print(f"  lowest: {min(usable):.4f}   highest: {max(usable):.4f}")
    else:
        print("\n  No cell type could be assessed, so no stability figures are available.")


if __name__ == "__main__":
    main()