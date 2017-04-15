#!/usr/bin/python2

import argparse
import os
import sys
import pandas as pd
import numpy as np


def create_output_file(outPath):
    '''
    holder
    '''
    output = '{}/yRNA_fragment_distribution.csv'.format(outPath)
    with open(output, 'w') as fo:
        fo.write('yRNA,loc,count,sample\n')
    return output


def load_data(fi):
    '''
    Load file
    '''
    return pd.read_csv(fi,
                       sep='\t')


def split_name_into_parts(df):
    '''
    Split the yRNA name into parts and subset df.
    '''
    df1 = df['Name'].str.split('_', expand = True)
    df1['yRNA'] = df1[0].str.split('-').str[1]
    df2 = pd.concat([df['Name'], df1['yRNA'], df1.iloc[ : , 1:-1], df['Count']], axis=1)
    df2.columns = ['Name', 'yRNA', 'start', 'end', 'length', 'strand', 'count']
    df2 = df2.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    return df2


def fragment_location_counts(df):
    '''
    holder
    '''
    frag_di = {}
    for index, row in df.iterrows():
        yRNA = row['yRNA']
        counts = row['count']

        if yRNA not in frag_di:
            frag_di[yRNA] = {l: 0 for l in range(0, row['length'] + 1)}

        for n in range(row['start'], row['end'] + 1):
            try:
                frag_di[yRNA][n] += counts
            except KeyError:
                pass
    return frag_di


def write_fragment_locations(frag_loc, output, name):
    '''
    holder
    '''
    with open(output, 'a+') as fo:
        for yRNA, locs in frag_loc.iteritems():
            for loc, count in locs.iteritems():
                fo.write('{},{},{},{}\n'.format(yRNA, loc, count, name))


def location_within_yRNA(df):
    '''
    Add whether rna fragment is in the 5', middle, or 3' end of yRNA.
    '''
    type_di = {1 : "5' half", 0 : "middle", -1 : "3' half"} 

    df['type'] = np.where(df['end'] < df['length'] / 2, 1, 
                          np.where(df['start'] > df['length'] / 2, -1, 0))
    df['type'] = np.where(df['strand'] == '-', df['type'] * -1,
                          np.where(df['strand'] == 'R-', df['type'] * -1, df['type'] * 1))

    df['type'].replace(type_di, inplace = True)

    summ_df = df[['yRNA','type','count']].groupby(['yRNA', 'type']).sum()
    summ_df.head()
    return df, summ_df


def write_output_files(df, summ_df, name):
    '''
    holder
    '''
    df.to_csv('tRNA_half.csv')
    summ_df.reset_index().set_index(['yRNA','type']).sortlevel(0).to_csv('tRNA_half_summary.csv')


def make_R_figures(r_dir, output):
    '''
    holder
    '''
    os.system('Rscript --vanilla {}/yrnaFragDistGraph.R {}'.format(r_dir, output))
    pass


def main(outPath, samples):
    output = create_output_file(outPath)
    for s in samples:
        name = s.split('/')[-3]
        df = load_data(s)
        df = split_name_into_parts(df)
        frag_loc = fragment_location_counts(df)
        write_fragment_locations(frag_loc, output, name)
    print 'Output saved as: {}'.format(output)
    make_R_figures(os.path.dirname(__file__), output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description='Analyzed the length distribution of trimmed reads')
    parser.add_argument(
             'outPath', 
             action='store', 
             help='Path to where the output file will be located')
    parser.add_argument(
             'samples', 
             action='store', 
             nargs='+',
             help='Path to where the sample output folders are located')
    arg = parser.parse_args()
    main(arg.outPath, arg.samples)
