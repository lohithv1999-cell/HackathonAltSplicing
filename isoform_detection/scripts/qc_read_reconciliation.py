#!/usr/bin/env python3
"""QC metric 3b: Read-count reconciliation (closing the IC:i:0 question).

The IC:i:0 profiling established that unassigned reads are overwhelmingly absent
from kallisto's pseudoalignment output. This script settles WHY they are absent, by
tracing the read count through the two steps between the input BAM and kallisto's
scored output:

    tagged BAM  --(collate + fastq)-->  input FASTQ  --(kallisto quant)-->  scored reads

For one cell type it reports:
  - primary reads in the tagged BAM            (the reads we started with)
  - reads in the input FASTQ                    (what actually reached kallisto)
  - reads kallisto scored (in pseudoalign BAM)  (what kallisto could place)

and the two differences between them. The interpretation:

  * If FASTQ count is close to the tagged-BAM primary count, essentially all reads
    reached kallisto -> the "never pseudo-aligned" reads genuinely matched no
    assembled transcript (the biologically interesting case).

  * If the FASTQ count is much smaller than the tagged-BAM count, reads were lost in
    the BAM->FASTQ conversion, so part of "never pseudo-aligned" is really "never
    reached kallisto" and must be qualified.

  * If the FASTQ count is close to the scored count, kallisto scored essentially
    everything it received; if it is much larger, kallisto itself dropped reads
    (e.g. too short for its k-mer length).

FASTQ read counting: a gzipped FASTQ has four lines per read, so the read count is
the line count divided by four. We read the gzip directly in Python so nothing needs
to be uncompressed on disk.
"""

import pysam
import argparse
import csv
import gzip


def count_primary_reads(bam_path):
    """Primary reads only (skip secondary/supplementary/unmapped)."""
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            n += 1
    return n


def count_all_records(bam_path):
    """Every record in a BAM (used for the pseudoalignment BAM, counting primaries)."""
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            n += 1
    return n


def count_fastq_reads(fastq_path):
    """Reads in a gzipped FASTQ = lines / 4."""
    lines = 0
    with gzip.open(fastq_path, "rt") as fh:
        for _ in fh:
            lines += 1
    return lines // 4


def main():
    parser = argparse.ArgumentParser(
        description="Reconcile read counts through BAM -> FASTQ -> kallisto to explain unassigned reads"
    )
    parser.add_argument("-t", "--tagged", required=True, help="Final IC-tagged BAM")
    parser.add_argument("-f", "--fastq", required=True, help="Input FASTQ (gzipped) fed to kallisto")
    parser.add_argument("-p", "--pseudobam", required=True, help="kallisto pseudoalignments.bam")
    parser.add_argument("-o", "--output", required=True, help="Output TSV report")
    parser.add_argument("-s", "--sample", required=True, help="Sample / cell type name")
    args = parser.parse_args()

    tagged_reads = count_primary_reads(args.tagged)
    fastq_reads = count_fastq_reads(args.fastq)
    scored_reads = count_all_records(args.pseudobam)

    # Where reads are lost, if anywhere
    lost_in_fastq_conversion = tagged_reads - fastq_reads
    lost_in_kallisto = fastq_reads - scored_reads

    def pct(x, base):
        return round(x / base * 100, 2) if base > 0 else 0.0

    with open(args.output, "w", newline="") as out_file:
        w = csv.writer(out_file, delimiter="\t")
        w.writerow([
            "Sample",
            "Tagged_BAM_reads",          # primary reads we started with
            "FASTQ_reads",               # reads that reached kallisto
            "Kallisto_scored_reads",     # reads kallisto placed (got a probability)
            "Lost_in_FASTQ_conversion",  # tagged - fastq
            "Lost_in_kallisto",          # fastq - scored
            "Pct_reached_kallisto",      # fastq / tagged
            "Pct_scored_of_reached",     # scored / fastq
        ])
        w.writerow([
            args.sample,
            tagged_reads,
            fastq_reads,
            scored_reads,
            lost_in_fastq_conversion,
            lost_in_kallisto,
            pct(fastq_reads, tagged_reads),
            pct(scored_reads, fastq_reads),
        ])

    print(f"[{args.sample}] tagged={tagged_reads:,}  fastq={fastq_reads:,}  "
          f"scored={scored_reads:,}  | lost_in_conversion={lost_in_fastq_conversion:,}  "
          f"lost_in_kallisto={lost_in_kallisto:,}")


if __name__ == "__main__":
    main()