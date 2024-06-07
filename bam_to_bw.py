import subprocess as sp
import logging
import pysam
import pyranges as pr
import argparse
import numpy as np
import os
import warnings

warnings.filterwarnings("ignore")

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script generates BigWig file from a BAM file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--bam_file", type=str, default=None)
    parser.add_argument("--peak_file", type=str, default=None)
    parser.add_argument("--extend_size", type=int, default=0)
    parser.add_argument("--forward_shift", type=int, default=4)
    parser.add_argument("--reverse_shift", type=int, default=-4)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    parser.add_argument("--chrom_size_file", type=str, default=None)

    return parser.parse_args()


def get_count(
    chrom: str = None,
    start: int = None,
    end: int = None,
    forward_shift: int = None,
    reverse_shift: int = None,
    extend_size: int = 0,
    bam: pysam.Samfile = None,
) -> np.array:
    """
    Get Tn5 cutting count from specific genomic region

    Parameters
    ----------
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    bam : pysam.Samfile
        BAM file
    """

    signal = np.zeros(shape=(end - start))

    for read in bam.fetch(reference=chrom, start=start, end=end):
        # cut counts
        if read.is_reverse:
            cut_site = read.reference_end + reverse_shift
        else:
            cut_site = read.reference_start + forward_shift

        if start <= cut_site < end:
            if extend_size > 0:
                _start = max(0, cut_site - start - extend_size)
                _end = min(cut_site - start + extend_size, end)
                signal[_start: _end] += 1
            else:
                signal[cut_site - start] += 1

    return signal


def main():
    args = parse_args()

    bam = pysam.Samfile(args.bam_file, "rb")

    logging.info(f"Loading genomic regions from {args.peak_file}")
    grs = pr.read_bed(args.peak_file)
    grs = grs.merge()

    logging.info(f"Total of {len(grs)} regions")

    wig_filename = os.path.join(args.out_dir, "{}.wig".format(args.out_name))
    bw_filename = os.path.join(args.out_dir, "{}.bw".format(args.out_name))

    # Open a new bigwig file for writing
    with open(wig_filename, "w") as f:
        for chrom, start, end in zip(grs.Chromosome, grs.Start, grs.End):
            signal = get_count(chrom=chrom, start=start, end=end,
                               forward_shift=args.forward_shift, 
                               reverse_shift=args.reverse_shift,
                               extend_size=args.extend_size,
                               bam=bam)

            f.write(f"fixedStep chrom={chrom} start={start+1} step=1\n")
            f.write("\n".join(str(e) for e in signal))
            f.write("\n")

    # convert to bigwig file
    logging.info("Conveting wig to bigwig!")
    sp.run(["wigToBigWig", wig_filename, args.chrom_size_file, bw_filename])
    os.remove(wig_filename)
    logging.info("Done!")


if __name__ == "__main__":
    main()
