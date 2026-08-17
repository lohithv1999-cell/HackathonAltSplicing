import pysam
import csv
import sys
import os

# This script takes a mixed BAM file and a CSV mapping barcodes to cell types, and splits the BAM file into separate BAM files for each cell type based on the barcode information in the BAM file.
def split_multiple_bams(bam_paths, csv_paths, output_dir, col_barcode="barcode", col_celltype="cell_type", col_bam_file="bam_file"):
    mapping = {}

    for csv_path in csv_paths:
        with open(csv_path, mode='r') as csv_file:
            reader = csv.DictReader(csv_file)
            has_bam_col = col_bam_file in reader.fieldnames

            for row in reader:
                cell_type = row[col_celltype].replace(" ", "_") # updated as per snakemake
                barcode = row[col_barcode]

            # If the user specified which BAM this barcode belongs to, we will use that information to create a more specific mapping. If not, we will assume that the barcode is unique across all BAM files.
                if has_bam_col:
                # If the CSV has a 'bam_file' column, we use it to create a mapping that includes both the BAM file name and the barcode. This allows us to handle multiple BAM files with potentially overlapping barcodes.
                    bam_name = os.path.basename(row[col_bam_file])
                    mapping[(bam_name, barcode)] = cell_type
                else:
                    mapping[barcode] = cell_type

    output_files = {}
    unique_cell_types = set(mapping.values())

    # Need to ensure that all BAM files have the same reference sequence dictionary. 
    # Then check the references and lengths of the first BAM file against all others.

    with pysam.AlignmentFile(bam_paths[0], "rb") as first_bam:
        header = first_bam.header.to_dict()
        reference = (tuple(first_bam.references), tuple(first_bam.lengths))

    seen_rgs = {rg["ID"] for rg in header.get("RG", [])}
    for bam_path in bam_paths[1:]:
        with pysam.AlignmentFile(bam_path, "rb") as extra_bam:
            if (tuple(extra_bam.references), tuple(extra_bam.lengths)) != reference:
                raise ValueError(
                    f"{bam_path} is aligned to a different reference from "
                    f"{bam_paths[0]}. All input BAMs must share a sequence dictionary."
                )
            for rg in extra_bam.header.to_dict().get("RG", []):
                if rg["ID"] in seen_rgs:
                    print(f"  warning: read group {rg['ID']} appears in more than one "
                          f"input; reads from {os.path.basename(bam_path)} will be "
                          f"indistinguishable from the earlier one")
                else:
                    header.setdefault("RG", []).append(rg)
                    seen_rgs.add(rg["ID"])

    print(f"Header carries {len(seen_rgs)} read group(s): {', '.join(sorted(seen_rgs))}")

    for ct in unique_cell_types:
        out_name = os.path.join(output_dir, f"{ct}.bam")
        output_files[ct] = pysam.AlignmentFile(out_name, "wb", header=header)

    # Iterate through each BAM file and process reads
    for bam_path in bam_paths:
        bam_name = os.path.basename(bam_path)
        print(f"Processing BAM file: {bam_name}")

        mixed_bam = pysam.AlignmentFile(bam_path, "rb")

    # Process each read in the BAM file
        for read in mixed_bam.fetch(until_eof=True):
            if read.has_tag('CB'):
                barcode = read.get_tag('CB')

            # Checking mapping dictionary
                cell_type = None
                if has_bam_col and (bam_name, barcode) in mapping:
                    cell_type = mapping[(bam_name, barcode)]
                elif not has_bam_col and barcode in mapping:
                    cell_type = mapping[barcode]

            # If a match is found, write it to cell specific file
                if cell_type:
                    output_files[cell_type].write(read)
        
        # closes the current BAM file before moving to next one
        mixed_bam.close()
                
# Close all output BAM files after processing
    for f in output_files.values():
        f.close()
    print("Finished splitting BAM files.")


if __name__ == "__main__":
# It reads the input BAM files and CSV mapping file from the Snakemake configuration, creates the output directory if it doesn't exist, and calls the function to split the BAM files.
    input_bams = snakemake.input.bams
    input_csv = snakemake.input.csvs
    output_dir = snakemake.params.output_dir

    p_barcode = snakemake.params.col_barcode
    p_cell_type = snakemake.params.col_celltype
    p_bam_file = snakemake.params.col_bam_file

    os.makedirs(output_dir, exist_ok=True)
    split_multiple_bams(input_bams, input_csv, output_dir, p_barcode, p_cell_type, p_bam_file)