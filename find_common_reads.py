import warnings

warnings.filterwarnings("ignore")

import os
import argparse
import pysam
import pandas as pd
import logging

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script find common reads in two bam files",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bam1", 
                        type=str, default=None, help="input BAM file")
    parser.add_argument("--bam2", 
                        type=str, default=None, help="input BAM file")
    
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
        help=(
            "If specified all output files will be written to that directory. \n"
            "Default: the current working directory"
        ),
    )
    parser.add_argument(
        "--out_name",
        type=str,
        default=None,
        help=("Names for output file. Default: output"),
    )

    return parser.parse_args()


def get_read_names(bam_path):
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    read_names = {read.query_name for read in bam_file.fetch()}
    bam_file.close()
    return read_names


def main():
    args = parse_args()

    read_names1 = get_read_names(args.bam1)
    read_names2 = get_read_names(args.bam2)
    
    common_names = list(read_names1.intersection(read_names2))
    
    df = pd.DataFrame(data={'reads': common_names})
    df.to_csv(f'{args.out_dir}/{args.name}.txt')