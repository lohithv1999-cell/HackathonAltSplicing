#!/usr/bin/env python3
"""QC metric 2: Parameter sensitivity (the ZW threshold sweep).

This script asks: if we changed the confidence cutoff, how many reads would still
be kept? It opens kallisto's pseudoalignment BAM (one cell type), reads the real
ZW probability that kallisto assigned to each read, and counts how many reads pass
each of several cutoffs.

If the pass rate barely changes as the cutoff moves, the pipeline is robust to the
threshold choice. If it changes a lot, the choice of cutoff matters and needs
justifying. Either way, this measures it directly from the real ZW values rather
than guessing.

A note on what is being counted. This file only contains reads that kallisto
actually pseudo-aligned, so the percentages here are out of "kallisto-aligned
reads", NOT out of the whole sequencing library. That is a different denominator
from the tag-completeness metric, which is why the two numbers are not expected to
match. We also skip secondary, supplementary, and unmapped records so each read is
counted once.
"""

import pysam
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Sweep kallisto ZW probabilities across several thresholds"
    )
    parser.add_argument("-b", "--bam", required=True, help="Input pseudoalignments.bam from kallisto")
    parser.add_argument("-o", "--output", required=True, help="Output TSV report")
    parser.add_argument("-s", "--sample", required=True, help="Sample / cell type name")
    args = parser.parse_args()

    # The cutoffs we want to test.
    thresholds = [0.50, 0.75, 0.80, 0.90, 0.95]

    # kallisto stores ZW as a 32-bit float, so a probability that is "really" 0.80
    # can come back as 0.7999999 and just miss a strict >= 0.80 test. We allow a
    # tiny tolerance so reads sitting exactly on a cutoff are counted as passing,
    # which is what a reader intuitively expects. With millions of reads this only
    # affects the handful sitting exactly on a boundary.
    epsilon = 1e-6

    # One running tally per cutoff, all starting at zero.
    passing = {t: 0 for t in thresholds}

    total_reads = 0     # reads that had a ZW probability we could read
    no_zw_reads = 0     # reads with no ZW tag (should be rare; tracked so it isn't hidden)

    with pysam.AlignmentFile(args.bam, "rb") as sam:
        for read in sam.fetch(until_eof=True):

            # Count each real read only once.
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            # We can only judge a read if it has a ZW probability.
            if not read.has_tag("ZW"):
                no_zw_reads += 1
                continue

            total_reads += 1
            zw_prob = read.get_tag("ZW")

            # For every cutoff, record whether this read would pass it.
            for t in thresholds:
                if zw_prob >= t - epsilon:
                    passing[t] += 1

    with open(args.output, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        writer.writerow([
            "Sample",
            "Threshold",
            "Passing_Reads",
            "Total_Reads",     # total reads that had a ZW value (the denominator)
            "Retention_Rate",  # percent of those reads that pass this cutoff
            "Reads_No_ZW",     # reads skipped because they had no ZW tag
        ])
        for t in thresholds:
            retention = (passing[t] / total_reads * 100) if total_reads > 0 else 0.0
            writer.writerow([
                args.sample,
                t,
                passing[t],
                total_reads,
                round(retention, 2),
                no_zw_reads,
            ])


if __name__ == "__main__":
    main()