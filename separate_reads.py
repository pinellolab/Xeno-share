import pandas as pd
import pysam
import argparse
import logging

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Saves reads below a alignment threshold and discards all others')
    parser.add_argument('--human_bam', type=str, default=None,
                        help="BAM file containing reads that mapped to human genome")
    parser.add_argument('--mouse_bam', type=str, default=None,
                        help="BAM file containing reads that mapped to mouse genome")
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--out_name', type=str, default=None)

    return parser.parse_args()


def count_reads(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as file:
        count = sum(1 for _ in file)
    return count


def main():
    args = parse_args()

    human_bam = pysam.AlignmentFile(args.human_bam, mode='rb')
    mouse_bam = pysam.AlignmentFile(args.mouse_bam, mode='rb')

    # check if they have same number of reads
    num_reads_human = sum(1 for _ in human_bam)
    num_reads_mouse = sum(1 for _ in mouse_bam)

    logging.info(f'Number of reads in human bam file: {num_reads_human}')
    logging.info(f'Number of reads in mouse bam file: {num_reads_mouse}')

    human_bam = pysam.AlignmentFile(args.human_bam, mode='rb')
    mouse_bam = pysam.AlignmentFile(args.mouse_bam, mode='rb')

    out_bam1 = pysam.AlignmentFile(f"{args.out_dir}/{args.out_name}_human.bam",
                                   "wb", template=human_bam)
    out_bam2 = pysam.AlignmentFile(f"{args.out_dir}/{args.out_name}_mouse.bam",
                                   "wb", template=mouse_bam)
    out_bam3 = pysam.AlignmentFile(f"{args.out_dir}/{args.out_name}_human_ambiguous.bam",
                                   "wb", template=human_bam)
    out_bam4 = pysam.AlignmentFile(f"{args.out_dir}/{args.out_name}_mouse_ambiguous.bam",
                                   "wb", template=mouse_bam)

    logging.info(f'Classifing reads')
    n_human, n_mouse, n_ambiguous = 0, 0, 0
    for human_read, mouse_read in zip(human_bam,
                                      mouse_bam):
        if human_read.qname != mouse_read.qname:
            logging.error('Reads have different name!')

        # compare alignment score
        try:
            human_as = human_read.get_tag('AS')
        except KeyError:
            human_as = -10000

        try:
            mouse_as = mouse_read.get_tag('AS')
        except KeyError:
            mouse_as = -10000

        if human_as > mouse_as:
            # write to human BAM file
            out_bam1.write(human_read)
            n_human += 1
        elif human_as < mouse_as:
            # write to mouse BAM file
            out_bam2.write(mouse_read)
            n_mouse += 1
        else:
            # write to ambiguous BAM file
            out_bam3.write(human_read)
            out_bam4.write(mouse_read)
            n_ambiguous += 1

    logging.info(f'Number of human reads: {n_human}')
    logging.info(f'Number of mouse reads: {n_mouse}')
    logging.info(f'Number of ambiguous reads: {n_ambiguous}')
    
    df = pd.DataFrame(data={'human': n_human, 
                            'mouse': n_mouse,
                            'ambiguous': n_ambiguous}, index=[0])
    df.to_csv(f'{args.out_dir}/{args.out_name}_summary.csv', 
              index=False)
    logging.info('Done!')
    

if __name__ == "__main__":
    main()
