#!/usr/bin/env python3
"""QC metric 3: IC:i:0 profiling (classification of unassigned reads).

This script answers a question the completeness and sensitivity metrics raised but
could not settle on their own: WHY is a read unassigned?

An unassigned read (IC:i:0 in the final BAM) can end up that way for two genuinely
different reasons, and until now we have only been able to infer the split
arithmetically. This script measures it directly, by comparing read identities
(QNAMEs) between the final tagged BAM and kallisto's pseudoalignment BAM.

For a chosen cell type it classifies every primary read in the final BAM into:

  - assigned            : IC tag > 0 (placed on a real isoform)
  - scored_but_low      : IC:i:0, AND the read IS present in kallisto's
                          pseudoalignment output -> kallisto scored it but the
                          probability fell below the assignment threshold. These
                          are the genuinely ambiguous reads.
  - not_pseudoaligned   : IC:i:0, AND the read is ABSENT from kallisto's
                          pseudoalignment output -> kallisto never scored it. These
                          are the reads that matched no assembled transcript well
                          enough to receive a probability. This is the group that
                          may point to splicing not captured by the current
                          assembly / annotation.

By reporting these three numbers per cell type, we replace the earlier inference
("of the order of seventeen million reads must be un-pseudoaligned") with a measured
fact, and we isolate the biologically interesting group for closer inspection.

Only primary reads are counted (secondary/supplementary/unmapped skipped), so each
read is counted once, consistent with the other QC scripts.

Memory note: the set of QNAMEs from the pseudoalignment BAM is held in memory. For
very deep cell types this can be large. The script therefore stores a hash of each
QNAME rather than the full string, which keeps memory manageable while remaining
effectively collision-free for this purpose. If you would rather store the exact
strings, set USE_HASH = False below.
"""

import pysam
import argparse
import csv

USE_HASH = True  # store hashed QNAMEs to save memory; set False to store exact strings


def load_pseudoaligned_ids(bam_path):
    """Return a set of read identities that appear in kallisto's pseudoalignment BAM.

    A read is considered 'pseudoaligned' if it appears here as a primary record.
    We key on QNAME (the read name), which is what links the same read across the
    two BAM files.
    """
    ids = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            name = read.query_name
            ids.add(hash(name) if USE_HASH else name)
    return ids


def main():
    parser = argparse.ArgumentParser(
        description="Classify unassigned (IC:i:0) reads into scored-but-low vs never-pseudoaligned"
    )
    parser.add_argument("-t", "--tagged", required=True,
                        help="Final IC-tagged BAM (the wp3_add_ic_tag output)")
    parser.add_argument("-p", "--pseudobam", required=True,
                        help="kallisto pseudoalignments.bam for the same cell type")
    parser.add_argument("-o", "--output", required=True, help="Output TSV report")
    parser.add_argument("-s", "--sample", required=True, help="Sample / cell type name")
    args = parser.parse_args()

    # Step 1: learn which reads kallisto actually scored.
    scored_ids = load_pseudoaligned_ids(args.pseudobam)

    # Step 2: walk the final tagged BAM and classify each primary read.
    assigned = 0
    scored_but_low = 0     # IC:i:0 but present in kallisto's output
    not_pseudoaligned = 0  # IC:i:0 and absent from kallisto's output
    missing_tag = 0        # no IC tag at all (expected to be zero for this pipeline)

    with pysam.AlignmentFile(args.tagged, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            try:
                ic_value = read.get_tag("IC")
            except KeyError:
                missing_tag += 1
                continue

            if ic_value > 0:
                assigned += 1
                continue

            # ic_value == 0: an unassigned read. Was it scored by kallisto or not?
            key = hash(read.query_name) if USE_HASH else read.query_name
            if key in scored_ids:
                scored_but_low += 1
            else:
                not_pseudoaligned += 1

    total = assigned + scored_but_low + not_pseudoaligned + missing_tag
    total_unassigned = scored_but_low + not_pseudoaligned + missing_tag

    def pct(x):
        return round(x / total * 100, 2) if total > 0 else 0.0

    with open(args.output, "w", newline="") as out_file:
        w = csv.writer(out_file, delimiter="\t")
        w.writerow([
            "Sample",
            "Total_Reads",
            "Assigned",
            "Unassigned_Total",
            "Scored_But_Low",        # IC:i:0, kallisto scored it (ambiguous)
            "Not_Pseudoaligned",     # IC:i:0, kallisto never scored it (off-assembly)
            "Missing_Tag",           # no IC tag at all
            "Pct_Assigned",          
            "Pct_Scored_But_Low",
            "Pct_Not_Pseudoaligned",
        ])
        w.writerow([
            args.sample,
            total,
            assigned,
            total_unassigned,
            scored_but_low,
            not_pseudoaligned,
            missing_tag,
            pct(assigned),
            pct(scored_but_low),
            pct(not_pseudoaligned),
        ])

    # A short human-readable summary to stdout as well.
    print(f"[{args.sample}] total={total:,}  assigned={assigned:,}  "
          f"scored_but_low={scored_but_low:,}  not_pseudoaligned={not_pseudoaligned:,}  "
          f"missing_tag={missing_tag:,}")


if __name__ == "__main__":
    main()