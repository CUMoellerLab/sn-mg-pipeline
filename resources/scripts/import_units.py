#!/usr/bin/env python

import click
import re
import pandas as pd
from glob import glob
from os.path import join, abspath


def parse_seq_dir(dir_fp,
                  extensions=['*.fastq',
                              '*.fastq.gz',
                              '*.fq',
                              '*.fq.gz']):

    seq_files = []
    for ext in extensions:
        seq_files.extend(glob(join(dir_fp, ext)))

    return(seq_files)


@click.command()
@click.option('-d', '--data_frame', 'df_path',
              required=True,
              help='path to data frame CSV')
@click.option('-m', '--match_column', 'match_column',
              required=True,
              help='name of column with search values')
@click.option('-n', '--name_column', 'name_column',
              required=False,
              help='name of column with sample names')
@click.option('-s', '--sequence_dir', 'seq_dirs',
              multiple=True,
              help='path to directory with sequences')
@click.option('-o', '--samples_fp', 'samples_fp',
              required=False,
              default='samples.txt',
              help='output samples.txt filepath')
@click.option('-u', '--units_fp', 'units_fp',
              required=False,
              default='units.txt',
              help='output units.txt filepath')
def import_df(df_path, match_column, seq_dirs, name_column,
              samples_fp, units_fp):
    """
    associates sample names from a dataframe column with f and r files
    """

    df = pd.read_csv(df_path)
    seq_fps = {}
    for seq_dir in seq_dirs:
        seq_fps[seq_dir] = parse_seq_dir(seq_dir)

    seq_tuples = []
    for i, row in df.iterrows():
        r1_pattern = re.compile('.+%s.+R1.+' % row[match_column])
        r2_pattern = re.compile('.+%s.+R2.+' % row[match_column])

        if name_column:
            name = row[name_column]
        else:
            name = i

        unit = 0
        for seq_dir in seq_fps:
            unit += 1

            unit_name = 'unit_%s' % unit
            r1 = None
            r2 = None
            for j, s in enumerate(seq_fps[seq_dir]):
                if re.match(r1_pattern, s):
                    r1 = abspath(s)
                if re.match(r2_pattern, s):
                    r2 = abspath(s)
                if r1 and r2:
                    seq_tuples.append((name, unit_name, r1, r2))
                    break

    with open(samples_fp, 'w') as f:
        f.write('Sample\n')
        for sample, _, _, _ in seq_tuples:
            f.write('%s\n' % sample)

    with open(units_fp, 'w') as f:
        f.write('Sample\tUnit\tR1\tR2\n')
        for sample, unit, r1, r2 in seq_tuples:
            f.write('%s\n' % '\t'.join([sample, unit, r1, r2]))


if __name__ == '__main__':
    import_df()
