#!/usr/bin/python2

import argparse
import datetime
import glob
import os
import sys
from yRNA_bin import yRNA_lenDist, \
                     yRNA_fragDist


def check_output_directory_path(outPath):
    '''
    holder
    '''
    if not os.path.isdir(outPath):
        print '\nERROR: The input directory does not exist'
        print 'Check to make sure path is correct\n'
        sys.exit()


def get_yRNA_summary_file_list(outPath):
    '''
    holder
    '''
    yRNAfiles = [y for y in glob.glob('{}/*/output/*yRNA.txt'.format(outPath))]
    if len(yRNAfiles) == 0:
        print '\nERROR: No TAB_3p_summary_yRNA.txt files detected'
        print 'Check to make sure output path correct and yRNA output exists'
        print 'in each samples output folder\n'
        sys.exit()
    return yRNAfiles


def create_output_folder(outPath):
    '''
    Create a uniquely named output folder for final results
    '''
    now = datetime.datetime.now()
    c = 1
    while True:
        d = '{}/{}_{}_{}_tRNAsummary_{}/'.format(outPath, now.year, now.month, now.day, c)
        try:
            os.makedirs(d)
            break
        except OSError:
            c += 1
    return d


def yRNA_length_distribution(outPath, samples):
    '''
    holder
    '''
    yRNA_lenDist.main(outPath, samples)


def yRNA_fragment_distribution(outPath, samples):
    '''
    holder
    '''
    yRNA_fragDist.main(outPath, samples)


def main(out_dir):
    check_output_directory_path(out_dir)
    yRNAfiles = get_yRNA_summary_file_list(out_dir)
    outPath = create_output_folder(out_dir)
    yRNA_length_distribution(outPath, yRNAfiles)
    yRNA_fragment_distribution(outPath, yRNAfiles)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description='Runs all the final summary scripts')
    parser.add_argument(
            'out_dir',
            action='store',
            help='Path to output directory of miRquant')
    arg = parser.parse_args()
    main(arg.out_dir)
