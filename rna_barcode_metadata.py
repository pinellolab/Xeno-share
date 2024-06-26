#!/usr/bin/env python3

"""
This script takes in a bam file, and outputs a txt file containing the number of
total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode.
"""

import argparse
import logging
import numpy as np
import pysam
from collections import defaultdict
from functools import partial


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get total reads, duplicate reads, UMIs, genes, and percent mitochondrial reads for each barcode from bam file")
    parser.add_argument("--bam_file", help="Filename for input bam file")
    parser.add_argument("--bai_file", help="Filename for bam index file")
    parser.add_argument("--genome", default="GRCh38")
    parser.add_argument("--barcode_metadata_file",
                        help="Filename for output barcode metadata txt file")
    parser.add_argument(
        "--subpool", help="Cellular subpool name", default=None, nargs="?")
    parser.add_argument(
        "--barcode_tag", help="BAM tag containing cell barcode", default="CB")

    return parser.parse_args()


def get_metrics(bam, barcode_tag="CB", subpool=None, genome=None):
    """
    Get barcode metrics from bam file; all counts are only for reads overlapping genes.
    Reported metrics are total counts, UMIs (one UMI counted per unique UMI-gene mapping),
    duplicate counts, genes, percent mitochondrial reads
    """

    # Reads per barcode is a vector of length 2.
    # [x_mito_reads, mito_reads]
    reads_per_barcode = defaultdict(partial(np.zeros, 2, dtype=int))
    genes_per_barcode = defaultdict(set)
    mito_genes_per_barcode = defaultdict(set)

    for read in bam:
        # skip read if not primary alignment (multimapper)
        if read.flag & 256:
            continue

        # get barcode; skip read if not present
        barcode = read.get_tag(barcode_tag)
        if barcode == "-":
            continue

        # get UMI; skip read if not present
        umi = read.get_tag("UB")
        if umi == "-":
            continue

        # get reads that mapped to genome;
        # if genome and genome not in read.reference_name:
        #     continue

        reads_per_barcode[barcode][0] += 1

        if "chrM" in read.reference_name:
            reads_per_barcode[barcode][1] += 1

        # get gene id; skip read if not present
        gene_id = read.get_tag("GX")
        if gene_id == "-":
            continue

        if "chrM" in read.reference_name:
            mito_genes_per_barcode[barcode].add(gene_id)
        else:
            genes_per_barcode[barcode].add(gene_id)

    # create list with barcodes and associated metrics
    barcode_metadata = []
    for barcode, reads_vector in reads_per_barcode.items():
        # Reminder that reads_vector is [total_reads, mito_reads].
        genes = len(genes_per_barcode[barcode])
        mito_genes = len(mito_genes_per_barcode[barcode])
        fraction_mitochondrial_reads = round(
            reads_vector[1]/reads_vector[0] * 100, 2)
        out_barcode = barcode + "_" + subpool if subpool else barcode

        metrics = list(map(str, [out_barcode,
                                 reads_vector[0],
                                 reads_vector[0]-reads_vector[1],
                                 reads_vector[1],
                                 genes+mito_genes,
                                 genes, mito_genes,
                                 fraction_mitochondrial_reads]))

        barcode_metadata.append(metrics)

    return barcode_metadata


def write_metadata_file(barcode_metadata, output_file):
    fields = ["barcode", "total_counts", "reads_non_mito",
              "reads_mito", "genes", "genes_non_mito",
              "genes_mito", "percent_mitochondrial"]

    with open(output_file, "w") as f:
        # write header
        f.write("\t".join(fields) + "\n")
        # write rows
        for metrics_list in barcode_metadata:
            f.write("\t".join(metrics_list) + "\n")


def main():
    # get arguments
    args = parse_arguments()
    bam_file = args.bam_file
    bai_file = args.bai_file
    subpool = args.subpool
    barcode_tag = args.barcode_tag
    barcode_metadata_file = args.barcode_metadata_file

    # load bam file
    # bam = pysam.AlignmentFile(bam_file, "rb", index_filename=bai_file)
    bam = pysam.AlignmentFile(bam_file, "rb")

    # get metrics for each barcode
    barcode_metadata = get_metrics(bam,
                                   barcode_tag,
                                   subpool,
                                   genome=args.genome)

    # write metadata file
    write_metadata_file(barcode_metadata, barcode_metadata_file)


if __name__ == "__main__":

    main()
