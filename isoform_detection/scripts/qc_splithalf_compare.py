#!/usr/bin/env python3
"""Compare the two halves of a split-half stability run.

Two things get measured. Assembly stability is how many isoforms turn up in both
halves, matched on intron chain with gffcompare rather than on name, since StringTie
numbers transcripts arbitrarily per run. Assignment stability is how well the read
counts agree across those matches.

Counts are compared raw, not as TPM, since TPM divides by effective length and that
depends on the fragment length setting.
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
    """Count transcripts, and return the ids of those with more than one exon.

    Single-exon transcripts have no intron chain, so gffcompare can't match them and
    they have to be kept out of the recovery denominator.
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
    """Pull exact-match pairs (class code "=") out of a gffcompare .tracking file.

    Columns: transfrag id, locus, reference gene|transcript, class code, then one
    column per input GTF. Locus comes back too, since transcripts sharing a locus
    compete for the same reads and it's worth being able to group by it.
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

            locus = fields[1]

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

            pairs.append((ref_tx, qry_tx, locus))
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


def transcript_coords(gtf_path):
    """Transcript coordinates from a GTF, for reporting unstable loci by position.

    gffcompare's XLOC ids are per-run, so they can't be compared across cell types.
    """
    coords = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if match:
                coords[match.group(1)] = (fields[0], int(fields[3]), int(fields[4]))
    return coords


def write_unstable_loci(path, sample, locus_a, locus_b, locus_txs, coords, top_n=20):
    """Write the loci whose counts disagree most between the halves.

    A handful of loci can account for most of the disagreement, and since they carry
    huge counts they drag the Pearson down while barely touching the rank
    correlation. Worth having them listed by coordinate rather than just a summary.
    """
    rows = []
    total = sum(abs(locus_a[k] - locus_b[k]) for k in locus_a)
    for locus in sorted(locus_a, key=lambda k: abs(locus_a[k] - locus_b[k]), reverse=True)[:top_n]:
        txs = locus_txs.get(locus, [])
        placed = [coords[t] for t in txs if t in coords]
        if placed:
            seqname = placed[0][0]
            start = min(p[1] for p in placed)
            end = max(p[2] for p in placed)
        else:
            seqname, start, end = "NA", "NA", "NA"
        diff = abs(locus_a[locus] - locus_b[locus])
        rows.append({
            "Sample": sample,
            "Seqname": seqname,
            "Start": start,
            "End": end,
            "N_Transcripts": len(txs),
            "Counts_A": round(locus_a[locus], 1),
            "Counts_B": round(locus_b[locus], 1),
            "Abs_Difference": round(diff, 1),
            "Pct_Of_Total_Discrepancy": round(100 * diff / total, 2) if total else 0.0,
        })

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows, total


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
    parser.add_argument("--unstable-loci", default=None,
                        help="Optional TSV listing the loci whose counts disagree most between halves")
    args = parser.parse_args()

    n_a, multi_a_ids = count_transcripts(args.gtf_a)
    n_b, multi_b_ids = count_transcripts(args.gtf_b)
    multi_a = len(multi_a_ids)
    multi_b = len(multi_b_ids)

    # Nothing to match on without intron chains, so say so instead of returning a
    # zero that looks like a result.
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
            "Loci_Compared": "NA",
            "Locus_Pearson": "NA",
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

    # Only count matches where both sides are multi-exon. gffcompare hands its "="
    # code to some single-exon transfrags too, and counting those against a
    # multi-exon denominator gave over 100%, which is how this got noticed.
    # The all-transcript figures are kept for context but understate things.
    multi_pairs = [(a, b, loc) for a, b, loc in pairs if a in multi_a_ids and b in multi_b_ids]
    n_matched_multi = len(multi_pairs)

    pct_a_recovered = 100 * n_matched / n_a if n_a else 0.0
    pct_b_recovered = 100 * n_matched / n_b if n_b else 0.0
    pct_a_multi_recovered = 100 * n_matched_multi / multi_a if multi_a else 0.0
    pct_b_multi_recovered = 100 * n_matched_multi / multi_b if multi_b else 0.0

    # Both transcript level and locus level, because they answer different things.
    # Where two isoforms overlap heavily the EM can give nearly all the reads to one
    # of them, and which one wins flips easily between halves. So a transcript can
    # swing from almost everything to almost nothing while its locus total barely
    # moves. Summing per locus first separates that from a real disagreement about
    # how much the gene was expressed.
    counts_a = read_abundance(args.abund_a)
    counts_b = read_abundance(args.abund_b)

    xs, ys = [], []
    locus_a, locus_b, locus_txs = {}, {}, {}
    for ref_tx, qry_tx, locus in pairs:
        if ref_tx in counts_a and qry_tx in counts_b:
            xs.append(counts_a[ref_tx])
            ys.append(counts_b[qry_tx])
            locus_a[locus] = locus_a.get(locus, 0.0) + counts_a[ref_tx]
            locus_b[locus] = locus_b.get(locus, 0.0) + counts_b[qry_tx]
            locus_txs.setdefault(locus, []).append(ref_tx)

    r_pearson = pearson(xs, ys)
    r_spearman = spearman(xs, ys)

    loci = sorted(locus_a)
    lxs = [locus_a[k] for k in loci]
    lys = [locus_b[k] for k in loci]
    r_locus = pearson(lxs, lys)

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
        "Loci_Compared": len(loci),
        "Locus_Pearson": round(r_locus, 4) if r_locus is not None else "NA",
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
    print(f"[{args.sample}] counts summed per locus: "
          f"Pearson={row['Locus_Pearson']} (n={len(loci):,} loci)")

    # Where the disagreement sits. If a few loci account for most of it, the Pearson
    # is being driven by outliers and the rank correlation is the fairer summary.
    if args.unstable_loci and locus_a:
        coords = transcript_coords(args.gtf_a)
        top_rows, total_diff = write_unstable_loci(
            args.unstable_loci, args.sample, locus_a, locus_b, locus_txs, coords)
        top_share = sum(r["Abs_Difference"] for r in top_rows)
        print(f"[{args.sample}] the {len(top_rows)} most discrepant loci account for "
              f"{100 * top_share / total_diff:.1f}% of all disagreement"
              if total_diff else f"[{args.sample}] no disagreement to report")


if __name__ == "__main__":
    main()