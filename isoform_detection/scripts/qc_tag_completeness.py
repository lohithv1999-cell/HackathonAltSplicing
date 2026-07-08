#!/usr/bin/env python3
"""QC metric 1: Tag completeness.

This script opens one final IC-tagged BAM file (one cell type) and asks a simple
question: of all the real reads in this file, how many actually received an
isoform-cluster (IC) tag?

It sorts every read into one of three groups:
  - assigned     : has an IC tag with a value greater than 0 (a real isoform)
  - unassigned   : has an IC tag equal to 0 (kallisto could not confidently place it)
  - missing_tag  : has no IC tag at all (the read never went through assignment)

Keeping "unassigned" and "missing_tag" separate matters: a read tagged IC:i:0 was
considered and judged ambiguous, whereas a read with no tag was never assessed.
They are different stories, so we count them separately.

We also count only PRIMARY reads. A single read can appear in a BAM several times
(secondary or supplementary alignments). If we counted all of those, one read would
be tallied many times and the totals would be wrong. We therefore skip secondary,
supplementary, and unmapped records and count each read once.
"""

import pysam
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Calculate IC tag completeness from a final tagged BAM"
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM file (tagged)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV report")
    parser.add_argument("-s", "--sample", required=True, help="Sample / cell type name")
    args = parser.parse_args()

    # Our three counters, all starting at zero.
    assigned = 0      # IC tag present and > 0
    unassigned = 0    # IC tag present and == 0
    missing_tag = 0   # no IC tag at all

    with pysam.AlignmentFile(args.input, "rb") as bam:
        for read in bam.fetch(until_eof=True):

            # Count each real read only once. Skip the extra alignment records
            # (secondary / supplementary) and anything that did not map.
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            # Try to read the IC tag. If it isn't there, this read was never tagged.
            try:
                ic_value = read.get_tag("IC")
            except KeyError:
                missing_tag += 1
                continue

            # The tag exists, so decide which bucket it belongs in.
            if ic_value > 0:
                assigned += 1
            else:
                unassigned += 1

    # "Unassigned" in the final report means everything that did not get a real
    # isoform: the explicit IC:i:0 reads plus the reads with no tag at all.
    # We still keep the breakdown so nothing is hidden.
    total = assigned + unassigned + missing_tag
    total_unassigned = unassigned + missing_tag

    percent_complete = (assigned / total * 100) if total > 0 else 0.0

    with open(args.output, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        # The first three columns match what the dashboard already expects.
        # The last two columns are the extra detail (the breakdown of "unassigned").
        writer.writerow([
            "Sample",
            "Total_Reads",
            "Assigned_Reads",
            "Unassigned_Reads",
            "Completeness_Percentage",
            "Unassigned_IC0",      # reads explicitly tagged IC:i:0
            "Unassigned_NoTag",    # reads with no IC tag at all
        ])
        writer.writerow([
            args.sample,
            total,
            assigned,
            total_unassigned,
            round(percent_complete, 2),
            unassigned,
            missing_tag,
        ])


if __name__ == "__main__":
    main()