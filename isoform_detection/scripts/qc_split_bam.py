#!/usr/bin/env python3
"""Split a BAM into two halves by hashing read names.

First step of the split-half stability check. Takes one pooled cell-type BAM and
deals every read into one of two output BAMs (half A or half B) based on a hash of
the read name. Each half then goes through the assembly and assignment chain
independently, and the two are compared.

The split is deterministic: the same read name always lands in the same half, so
re-running on the same input gives the same two halves. That makes the QC itself
repeatable, even though the halves deliberately contain different reads and so are
never expected to agree exactly.

Reads are single-end here, so there are no mates to keep together. Hashing the name
also means any secondary/supplementary records, which share a primary's read name,
travel to the same half.
"""

import argparse
import hashlib
import pysam


def half_for(read_name):
    """Return 0 or 1 for a read name, deterministically.

    md5 is used only as a stable, uniform hash. Its output does not depend on the
    Python process, unlike the built-in hash(), which is randomised per run. That
    process-independence is what makes the split reproducible across runs.
    """
    digest = hashlib.md5(read_name.encode()).digest()
    return digest[0] & 1


def main():
    parser = argparse.ArgumentParser(description="Split a BAM into two halves by read-name hash")
    parser.add_argument("-i", "--input", required=True, help="Input BAM (a pooled cell-type BAM)")
    parser.add_argument("-a", "--out_a", required=True, help="Output BAM for half A")
    parser.add_argument("-b", "--out_b", required=True, help="Output BAM for half B")
    args = parser.parse_args()

    n_a = 0
    n_b = 0

    with pysam.AlignmentFile(args.input, "rb") as infile:
        with pysam.AlignmentFile(args.out_a, "wb", template=infile) as out_a, \
             pysam.AlignmentFile(args.out_b, "wb", template=infile) as out_b:
            for read in infile:
                if half_for(read.query_name) == 0:
                    out_a.write(read)
                    n_a += 1
                else:
                    out_b.write(read)
                    n_b += 1

    total = n_a + n_b
    print(f"[split] {args.input}")
    print(f"[split] half A: {n_a:,}  half B: {n_b:,}  total: {total:,}")
    if total > 0:
        pct_a = 100 * n_a / total
        print(f"[split] half A is {pct_a:.2f}% of reads (expected ~50%)")


if __name__ == "__main__":
    main()