# Output column reference

Every quality-control file the workflow writes is tab-separated. Most carry one row per
cell type; where that is not so, it is said below. Columns whose meaning follows from
the name are left out, and the rest are explained.

Two conventions hold throughout. `Sample` is the cell type, or the sample name when the
workflow runs in bulk mode. Percentages are numbers between 0 and 100, not fractions.

---

## `qc/completeness/<cell>_completeness.tsv`

How much of a cell type's library was assigned an isoform.

| Column | Meaning |
|---|---|
| `Total_Reads` | Primary alignments in the tagged BAM. Secondary and supplementary alignments are left out, so this counts reads rather than alignment records. |
| `Assigned_Reads` | Reads whose `IC` tag is greater than zero, meaning they were attributed to a particular assembled transcript. |
| `Unassigned_Reads` | Reads not attributed to a transcript, for whatever reason. |
| `Completeness_Percentage` | `Assigned_Reads` as a percentage of `Total_Reads`. |
| `Unassigned_IC0` | Unassigned reads explicitly tagged `IC:i:0`. |
| `Unassigned_NoTag` | Unassigned reads carrying no `IC` tag at all. This should be zero, since the workflow tags every read; anything else means reads slipped past the tagging step. |

## `qc/sensitivity/<cell>_sensitivity.tsv`

How the assignment threshold affects what is kept. One row per threshold tested, so
several rows per cell type.

| Column | Meaning |
|---|---|
| `Threshold` | The `ZW` probability being tested. `ZW` is kallisto's estimate that a read has been assigned to the right transcript. |
| `Passing_Reads` | Reads whose `ZW` met or exceeded the threshold. |
| `Total_Reads` | Reads carrying a `ZW` value. Note that the denominator here is scored reads, not all reads. |
| `Retention_Rate` | `Passing_Reads` as a percentage of `Total_Reads`. This is retention among reads kallisto managed to score; reads it never pseudo-aligned have no `ZW` and do not appear here at all, so the figure is not comparable with completeness. |
| `Reads_No_ZW` | Reads in the pseudoalignment output that carry no `ZW` value. |

## `qc/ic0/<cell>_ic0.tsv`

Why each unassigned read failed. Two causes are distinguished, and on a correct run they
account for every unassigned read between them.

| Column | Meaning |
|---|---|
| `Assigned`, `Unassigned_Total` | The same quantities as `Assigned_Reads` and `Unassigned_Reads` in the completeness file, repeated so this file can be read on its own. |
| `Scored_But_Low` | Reads kallisto pseudo-aligned and gave a `ZW` probability that fell below the threshold. These reads matched transcripts, and the workflow was not confident enough to choose between them. |
| `Not_Pseudoaligned` | Reads kallisto could not place against any assembled transcript, so no probability exists. These reads matched nothing. |
| `Missing_Tag` | Reads carrying no `IC` tag. Expected to be zero, and a non-zero value points at the tagging step rather than at the data. |
| `Pct_Assigned` | A percentage of **all** reads. |
| `Pct_Scored_But_Low`, `Pct_Not_Pseudoaligned` | Percentages of **`Unassigned_Total`**, not of all reads, and the two sum to 100. The three percentage columns therefore do not share a denominator, which is worth keeping in mind when reading them side by side. |

## `qc/splithalf/<cell>_stability.tsv` and `ALL_stability.tsv`

Agreement between two independently processed halves of the same cell type. One file per
cell type, plus a collated `ALL_stability.tsv` carrying every row. Written only when the
split-half check is enabled.

Transcripts are matched between the halves on **intron chain** using gffcompare, not on
name. StringTie numbers transcripts arbitrarily on each run, so an identifier from one
assembly says nothing about the other. Only exact intron-chain matches count.

| Column | Meaning |
|---|---|
| `Transcripts_A`, `Transcripts_B` | Transcripts assembled in each half. |
| `MultiExon_A`, `MultiExon_B` | How many of those have more than one exon. A single-exon transcript has no intron chain, so it cannot be matched this way. |
| `Matched_Transcripts` | Transcript pairs that gffcompare matched exactly. |
| `Matched_MultiExon` | How many of those matches are multi-exon on both sides. |
| `Pct_MultiExon_A_Recovered` | `Matched_MultiExon` as a percentage of `MultiExon_A`. **This is the assembly-stability figure to read.** Numerator and denominator are both restricted to multi-exon transcripts, so they describe the same population. |
| `Pct_A_Recovered_In_B` | `Matched_Transcripts` as a percentage of `Transcripts_A`. Kept for context, though it understates stability: its denominator includes single-exon transcripts that could never have been matched. |
| `Pct_MultiExon_B_Recovered`, `Pct_B_Recovered_In_A` | The same two figures measured the other way round. In practice they differ from their counterparts by a fraction of a percentage point, since neither half is systematically richer than the other. |
| `Compared_Transcripts` | Matched pairs for which a read count exists in both halves, which is the number of points behind the correlations. |
| `Counts_Pearson`, `Counts_Spearman` | Correlation of kallisto's estimated counts across matched transcripts. Raw counts are used rather than TPM, because TPM divides by effective length and would tie the figure to the configured fragment length. Where the two correlations diverge, a few very abundant transcripts are pulling the Pearson value about, and the rank correlation is the fairer summary. |
| `Loci_Compared` | Loci represented among the matched transcripts. A locus here is a gffcompare grouping of overlapping transcripts, which usually corresponds to a gene but is not guaranteed to. |
| `Locus_Pearson` | The same correlation after summing matched transcripts within each locus. Set against `Counts_Pearson`, it separates disagreement about how much a gene was expressed from a reallocation of reads between that gene's isoforms. |
| `Total_Counts_A`, `Total_Counts_B` | Estimated counts summed over all transcripts in each half. |
| `Status` | Either `ok`, or `no_multi_exon_transcripts` where the assembly held nothing with an intron chain to match on. The latter usually means `library_type` has not been set correctly. No other values are emitted. |

## `qc/splithalf/unstable_loci/<cell>_unstable_loci.tsv`

The loci whose read counts disagree most between the halves, worst first. Written
alongside the stability table.

| Column | Meaning |
|---|---|
| `Seqname`, `Start`, `End` | Coordinates of the locus, spanning all its matched transcripts. Coordinates are given rather than gffcompare's locus identifiers, which are assigned afresh on each run and cannot be compared between cell types. |
| `N_Transcripts` | Matched transcripts in this locus. |
| `Counts_A`, `Counts_B` | Estimated counts summed across those transcripts in each half. |
| `Abs_Difference` | The absolute difference between the two. |
| `Pct_Of_Total_Discrepancy` | This locus's share of the summed absolute difference across every locus. A handful of loci accounting for a large share means the correlation is being driven by outliers rather than by noise spread thinly. |

## `qc/summary_statistics.tsv`

Cohort-level figures, written as `Metric` and `Value` pairs rather than one row per cell
type.

| Metric | Meaning |
|---|---|
| `n_cell_types` | Cell types included in the run. |
| `total_reads`, `total_assigned` | Summed across every cell type. |
| `pooled_assignment_rate_pct` | Assigned as a percentage of total across the pooled cohort. This is not the mean of the per-cell-type rates, which would give a cell type with a hundred thousand reads the same weight as one with a hundred million. |
| `completeness_corr_log_reads` | Pearson correlation between completeness and the base-ten logarithm of the read count per cell type. These are pooled reads per cell type rather than a per-cell sequencing depth: a cell type with many cells at low coverage and one with few cells at high coverage can carry the same total. |
| `completeness_min_pct`, `completeness_max_pct` | The extremes across cell types, each naming the cell type concerned. |
| `mean_retention_at_X_pct` | Mean retention at threshold X, averaged across cell types. Written once per threshold tested. |
| `retention_swing_pct` | Retention at the lowest threshold minus retention at the highest. A large value means the threshold materially changes what is kept, and a small one that it does not. |