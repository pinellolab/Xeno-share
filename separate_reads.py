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


def separate_atac(human_bam_file, mouse_bam_file, out_dir, out_name):
    human_bam = pysam.AlignmentFile(human_bam_file, mode='rb')
    mouse_bam = pysam.AlignmentFile(mouse_bam_file, mode='rb')

    out_bam1 = pysam.AlignmentFile(f"{out_dir}/{out_name}_human.bam",
                                   "wb", template=human_bam)
    out_bam2 = pysam.AlignmentFile(f"{out_dir}/{out_name}_mouse.bam",
                                   "wb", template=mouse_bam)
    out_bam3 = pysam.AlignmentFile(f"{out_dir}/{out_name}_human_ambiguous.bam",
                                   "wb", template=human_bam)
    out_bam4 = pysam.AlignmentFile(f"{out_dir}/{out_name}_mouse_ambiguous.bam",
                                   "wb", template=mouse_bam)

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

    return n_human, n_mouse, n_ambiguous


def separate(human_bam_file, mouse_bam_file, out_dir, out_name):
    human_bam = pysam.AlignmentFile(human_bam_file, mode='rb')
    mouse_bam = pysam.AlignmentFile(mouse_bam_file, mode='rb')

    human_reads, mouse_reads, ambiguous_reads = [], [], []
    
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
            human_reads.append(human_read.qname)
        elif human_as < mouse_as:
            mouse_reads.append(mouse_read.qname)
        else:
            ambiguous_reads.append(human_read.qname)

    df_human = pd.DataFrame({'read_name': human_reads})
    df_mouse = pd.DataFrame({'read_name': mouse_reads})
    df_ambiguous = pd.DataFrame({'read_name': ambiguous_reads})

    df_human.to_csv(f'{out_dir}/{out_name}_human.csv',
                    index=False, header=False, sep='\t')
    df_mouse.to_csv(f'{out_dir}/{out_name}_mouse.csv',
                    index=False, header=False, sep='\t')
    df_ambiguous.to_csv(
        f'{out_dir}/{out_name}_ambiguous.csv',
        index=False, header=False, sep='\t')

    return len(human_reads), len(mouse_reads), len(ambiguous_reads)


def main():
    args = parse_args()

    human_bam = pysam.AlignmentFile(args.human_bam, mode='rb')
    mouse_bam = pysam.AlignmentFile(args.mouse_bam, mode='rb')

    # check if they have same number of reads
    num_reads_human = sum(1 for _ in human_bam)
    num_reads_mouse = sum(1 for _ in mouse_bam)

    logging.info(f'Number of reads in human bam file: {num_reads_human}')
    logging.info(f'Number of reads in mouse bam file: {num_reads_mouse}')

    logging.info(f'Classifing reads')
    n_human, n_mouse, n_ambiguous = separate(human_bam_file=args.human_bam,
                                             mouse_bam_file=args.mouse_bam,
                                             out_dir=args.out_dir,
                                             out_name=args.out_name)

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
