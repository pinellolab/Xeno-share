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
    parser.add_argument('--human_reads', type=str, default=None)
    parser.add_argument('--mouse_reads', type=str, default=None)
    parser.add_argument('--ambiguous_reads', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--out_name', type=str, default=None)

    return parser.parse_args()


def main():
    args = parse_args()

    df_human = pd.read_csv(args.human_reads, sep="\t", header=None)
    df_mouse = pd.read_csv(args.mouse_reads, sep="\t", header=None)
    df_ambiguous = pd.read_csv(args.ambiguous_reads, sep="\t", header=None)

    df_human.columns = ['read']
    df_mouse.columns = ['read']
    df_ambiguous.columns = ['read']

    barcodes, species = [], []

    for read in df_human['read']:
        names = read.split('_')[1].split(',')
        barcode = names[0] + names[1] + names[2]

        barcodes.append(barcode)
        species.append('human')

    for read in df_mouse['read']:
        names = read.split('_')[1].split(',')
        barcode = names[0] + names[1] + names[2]

        barcodes.append(barcode)
        species.append('mouse')

    for read in df_ambiguous['read']:
        names = read.split('_')[1].split(',')
        barcode = names[0] + names[1] + names[2]

        barcodes.append(barcode)
        species.append('ambiguous')

    df_barcode = pd.DataFrame({'barcode': barcodes, 'species': species})
    df_count = df_barcode.groupby(
        ['barcode', 'species']).size().reset_index(name='count')
    df_count = df_count.pivot(index='barcode', columns='species',
                              values='count').fillna(0)

    df_count = df_count.astype({'ambiguous': 'int', 'human': 'int', 'mouse': 'int'})
    
    df_count.to_csv(f'{args.out_dir}/{args.out_name}.csv')
    logging.info('Done!')


if __name__ == "__main__":
    main()
