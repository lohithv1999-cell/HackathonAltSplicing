#!/usr/bin/env python3
"""Compare the two halves of a split-half stability run.

The split-half check runs each half of a cell type's reads through assembly and
pseudoalignment independently. This script compares the two results and reports how
far they agree.

Two things are measured, and they test different parts of the pipeline:

  Assembly stability   Did the two halves find the same isoforms? StringTie names
                       transcripts arbitrarily per run, so half A's STRG.1.1 and
                       half B's STRG.1.1 are unrelated. gffcompare matches them on
                       intron chain instead, and only exact matches (class code "=")
                       are counted.

  Assignment stability Across the isoforms found in both halves, do they receive
                       similar read counts? Reported as Pearson and Spearman
                       correlation of kallisto's est_counts.

Counts are compared raw rather than as TPM. TPM divides by effective length, which
depends on the configured fragment length, and comparing raw counts keeps this
measure independent of that parameter. The two halves carry near-identical read
totals in any case, so there is little to normalise away.

Note that the halves are expected to differ. They contain different reads by design,
so this measures stability, not reproducibility.
"""

import argparse
import csv
import math
import os
import re
import subprocess
import sys


def read_abundance(path):
    """Return {target_id: est_counts} from a kallisto abundance.tsv."""
    counts = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            counts[row["target_id"]] = float(row["est_counts"])
    return counts


def count_transcripts(gtf_path):
    """Count transcripts in a GTF and identify which have more than one exon.

    Returns (n_transcripts, set_of_multi_exon_transcript_ids).

    The multi-exon set matters because the whole comparison rests on it.
    gffcompare's exact-match class code is primarily an intron-chain match, so a
    single-exon transcript has no chain to match on. An assembly with no
    multi-exon transcripts cannot be compared at all, and the caller needs to be
    told that rather than handed a zero. Returning the identifiers rather than
    just a count also lets the matched pairs be filtered to multi-exon transcripts,
    so that the numerator and denominator describe the same population.

    Transcripts are counted from transcript features. Exons are tallied per
    transcript_id, since a transcript line does not itself say how many exons follow.
    """
    n_tx = 0
    exons_per_tx = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature == "transcript":
                n_tx += 1
            elif feature == "exon":
                match = re.search(r'transcript_id "([^"]+)"', fields[8])
                if match:
                    tx_id = match.group(1)
                    exons_per_tx[tx_id] = exons_per_tx.get(tx_id, 0) + 1

    multi_exon = {tx for tx, count in exons_per_tx.items() if count > 1}
    return n_tx, multi_exon


def run_gffcompare(gtf_a, gtf_b, prefix):
    """Run gffcompare with half A as reference and half B as query.

    Returns the path to the .tracking file, which lists each query transcript
    alongside the reference transcript it matched and the class code describing
    the relationship.
    """
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
    cmd = ["gffcompare", "-r", gtf_a, "-o", prefix, gtf_b]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"gffcompare failed for {gtf_a} vs {gtf_b}")
    return prefix + ".tracking"


def parse_tracking(tracking_path):
    """Pull exact-match transcript pairs out of a gffcompare .tracking file.

    Columns are: query transfrag id, locus id, reference gene|transcript, class
    code, then one column per input GTF holding that GTF's transcript for this
    transfrag. We keep only rows with class code "=", meaning the intron chains
    match exactly, and pull the transcript id from each side.
    """
    pairs = []
    with open(tracking_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            class_code = fields[3]
            if class_code != "=":
                continue

            # Reference column looks like "gene_id|transcript_id".
            ref_field = fields[2]
            if "|" not in ref_field:
                continue
            ref_tx = ref_field.split("|")[1]

            # Query column looks like "q1:gene|transcript|<numbers...>".
            qry_field = fields[4]
            if qry_field == "-":
                continue
            parts = qry_field.split("|")
            if len(parts) < 2:
                continue
            qry_tx = parts[1]

            pairs.append((ref_tx, qry_tx))
    return pairs


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


def rank(values):
    """Fractional ranks, averaging ties."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = mean_rank
        i = j + 1
    return ranks


def spearman(xs, ys):
    if len(xs) < 2:
        return None
    return pearson(rank(xs), rank(ys))


def write_row(path, row):
    """Write a single-row TSV. Used by both the normal and the not-applicable paths,
    so that a run which cannot be assessed still produces a file the collate step can
    read, carrying the reason in its Status column."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Compare two halves of a split-half stability run")
    parser.add_argument("--gtf-a", required=True, help="StringTie GTF for half A")
    parser.add_argument("--gtf-b", required=True, help="StringTie GTF for half B")
    parser.add_argument("--abund-a", required=True, help="kallisto abundance.tsv for half A")
    parser.add_argument("--abund-b", required=True, help="kallisto abundance.tsv for half B")
    parser.add_argument("--gffcompare-prefix", required=True, help="Output prefix for gffcompare files")
    parser.add_argument("--sample", required=True, help="Cell type name, for the output row")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    args = parser.parse_args()

    n_a, multi_a_ids = count_transcripts(args.gtf_a)
    n_b, multi_b_ids = count_transcripts(args.gtf_b)
    multi_a = len(multi_a_ids)
    multi_b = len(multi_b_ids)

    # gffcompare matches on intron chain, so an assembly with no multi-exon
    # transcripts offers nothing to match and would silently return zero. Say so
    # plainly instead, and stop before running a comparison that cannot work.
    if multi_a == 0 or multi_b == 0:
        write_row(args.output, {
            "Sample": args.sample,
            "Transcripts_A": n_a,
            "Transcripts_B": n_b,
            "MultiExon_A": multi_a,
            "MultiExon_B": multi_b,
            "Matched_Transcripts": "NA",
            "Matched_MultiExon": "NA",
            "Pct_MultiExon_A_Recovered": "NA",
            "Pct_MultiExon_B_Recovered": "NA",
            "Pct_A_Recovered_In_B": "NA",
            "Pct_B_Recovered_In_A": "NA",
            "Compared_Transcripts": "NA",
            "Counts_Pearson": "NA",
            "Counts_Spearman": "NA",
            "Total_Counts_A": "NA",
            "Total_Counts_B": "NA",
            "Status": "no_multi_exon_transcripts",
        })
        sys.stderr.write(
            f"[{args.sample}] Split-half stability cannot be assessed.\n"
            f"[{args.sample}] Half A has {multi_a:,} multi-exon transcripts of {n_a:,}; "
            f"half B has {multi_b:,} of {n_b:,}.\n"
            f"[{args.sample}] Isoforms are matched between halves on intron chain, so an\n"
            f"[{args.sample}] assembly without multi-exon transcripts cannot be compared.\n"
            f"[{args.sample}] This usually means StringTie ran without a reference annotation,\n"
            f"[{args.sample}] or the reads carry too little junction-spanning evidence to\n"
            f"[{args.sample}] assemble spliced transcripts. Supply an annotation and re-run.\n"
        )
        return

    tracking = run_gffcompare(args.gtf_a, args.gtf_b, args.gffcompare_prefix)
    pairs = parse_tracking(tracking)
    n_matched = len(pairs)

    # Assembly stability: how much of each half's assembly the other half recovered.
    #
    # The multi-exon figure is the meaningful one, and it is computed from matches
    # where both transcripts are multi-exon. Restricting the numerator this way
    # matters: gffcompare also assigns its exact-match code to some single-exon
    # transfrags, so counting every match against a multi-exon denominator can
    # exceed 100 per cent, which is how this was first noticed. Numerator and
    # denominator now describe the same population.
    #
    # The all-transcript figures are kept for context, but they understate stability,
    # since single-exon transcripts are largely unmatchable by intron chain and their
    # number varies with how fragmented the assembly happens to be.
    multi_pairs = [(a, b) for a, b in pairs if a in multi_a_ids and b in multi_b_ids]
    n_matched_multi = len(multi_pairs)

    pct_a_recovered = 100 * n_matched / n_a if n_a else 0.0
    pct_b_recovered = 100 * n_matched / n_b if n_b else 0.0
    pct_a_multi_recovered = 100 * n_matched_multi / multi_a if multi_a else 0.0
    pct_b_multi_recovered = 100 * n_matched_multi / multi_b if multi_b else 0.0

    # Assignment stability: correlate read counts across the matched isoforms only.
    counts_a = read_abundance(args.abund_a)
    counts_b = read_abundance(args.abund_b)

    xs, ys = [], []
    for ref_tx, qry_tx in pairs:
        if ref_tx in counts_a and qry_tx in counts_b:
            xs.append(counts_a[ref_tx])
            ys.append(counts_b[qry_tx])

    r_pearson = pearson(xs, ys)
    r_spearman = spearman(xs, ys)

    total_a = sum(counts_a.values())
    total_b = sum(counts_b.values())

    row = {
        "Sample": args.sample,
        "Transcripts_A": n_a,
        "Transcripts_B": n_b,
        "MultiExon_A": multi_a,
        "MultiExon_B": multi_b,
        "Matched_Transcripts": n_matched,
        "Matched_MultiExon": n_matched_multi,
        "Pct_MultiExon_A_Recovered": round(pct_a_multi_recovered, 2),
        "Pct_MultiExon_B_Recovered": round(pct_b_multi_recovered, 2),
        "Pct_A_Recovered_In_B": round(pct_a_recovered, 2),
        "Pct_B_Recovered_In_A": round(pct_b_recovered, 2),
        "Compared_Transcripts": len(xs),
        "Counts_Pearson": round(r_pearson, 4) if r_pearson is not None else "NA",
        "Counts_Spearman": round(r_spearman, 4) if r_spearman is not None else "NA",
        "Total_Counts_A": round(total_a, 1),
        "Total_Counts_B": round(total_b, 1),
        "Status": "ok",
    }
    write_row(args.output, row)

    print(f"[{args.sample}] transcripts A={n_a:,} B={n_b:,} "
          f"(multi-exon {multi_a:,} / {multi_b:,})")
    print(f"[{args.sample}] matched={n_matched:,} of which {n_matched_multi:,} multi-exon")
    print(f"[{args.sample}] assembly stability: {pct_a_multi_recovered:.1f}% of A's "
          f"multi-exon transcripts recovered in B ({pct_b_multi_recovered:.1f}% the other way)")
    print(f"[{args.sample}] counts on matched isoforms: "
          f"Pearson={row['Counts_Pearson']} Spearman={row['Counts_Spearman']} "
          f"(n={len(xs):,})")


if __name__ == "__main__":
    main()