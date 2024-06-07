import pandas as pd
import argparse
import logging
from Bio import SeqIO
import gzip

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_fastq', type=str, default=None)
    parser.add_argument('--human_reads', type=str, default=None)
    parser.add_argument('--mouse_reads', type=str, default=None)
    parser.add_argument('--human_fastq', type=str, default=None)
    parser.add_argument('--mouse_fastq', type=str, default=None)    
    return parser.parse_args()


def main():
    args = parse_args()

    df_human = pd.read_csv(args.human_reads, header=None, sep="\t")
    df_mouse = pd.read_csv(args.mouse_reads, header=None, sep="\t")
    
    df_human.columns = ['human_reads']
    df_mouse.columns = ['moues_reads']
    
    human_reads = set(df_human['human_reads'].tolist())
    mouse_reads = set(df_mouse['moues_reads'].tolist())

    in_f = gzip.open(args.in_fastq, "rt")
    out_f_human = gzip.open(args.human_fastq, "wt")
    out_f_mouse = gzip.open(args.mouse_fastq, "wt")

    logging.info('Spliting FASTQ file!')
    for record in SeqIO.parse(in_f, "fastq"):
        if record.id in human_reads:
            out_f_human.write(record.format("fastq"))
        elif record.id in mouse_reads:
            out_f_mouse.write(record.format("fastq"))

    in_f.close()
    out_f_human.close()
    out_f_mouse.close()

    logging.info('Done!')


if __name__ == "__main__":
    main()
