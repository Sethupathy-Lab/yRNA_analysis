#!/usr/bin/python2

usage='''
 Usage: python lenDist.py /path/to/sample_output_directory file_list

 Output saved as lenDist_yRNA.csv and lenDist_yRNA.png

 Description:
   Calculates the read length distribution post-trimming across all samples.
   Create a length distribultion graph for each of the samples.

'''
 
import argparse
import os
import sys


def calculate_yRNA_proportions(yRNA_di):
    '''
    holder
    '''
    tot = sum([sum(l.values()) for l in yRNA_di.values()])
    tot_count = {y: sum(l.values()) / tot * 100 for y, l in yRNA_di.iteritems()}
    return tot_count


def calculate_length_percentage(yRNA_di):
    '''
    Calculate the percentage of total reads for each
    read length on a yRNA-by-yRNA basis
    '''
    for yRNA, lengths in yRNA_di.iteritems():
        tot = sum(lengths.values())
        for len, count in lengths.iteritems():
            yRNA_di[yRNA][len] = count / tot * 100
    return yRNA_di


def read_lengths_dict(file):
    '''
    Load the length distribution data into a dictionary
    If the read length is greater than 60% of the yRNA length,
    don't include in the analysis.
    '''
    yRNA_di = {}
    lengths = {}
    yRNAs = {}
    with open(file, 'r') as fi:
        fi.next()
        for l in fi:
            l = l.split('\t')
            counts = float(l[5])
            name_li = l[0].split('_')
            yRNA = name_li[0].split('-')[-1]
            length = int(name_li[2]) - int(name_li[1])
            full_length = int(name_li[3])
           
            # Dont include if nearly a full yRNA
            if length > 50:
                print l
                continue

            if yRNA not in yRNA_di:
                yRNA_di[yRNA] = {}

            yRNAs[yRNA] = 1
            try:
                yRNA_di[yRNA][length] += counts
            except KeyError:
                yRNA_di[yRNA][length] = counts
            lengths[length] = 1

        tot_count = calculate_yRNA_proportions(yRNA_di)
        yRNA_di = calculate_length_percentage(yRNA_di)
    return yRNA_di, lengths, yRNAs, tot_count


def process_yRNA_files(samples):
    '''
    For each TAB_3p_summary_yRNA file, get a summary of the
    read length distribution
    '''
    out_di = {f.split('/')[-3]: {} for f in samples}
    tot_counts = {f.split('/')[-3]: {} for f in samples}
    lengths = {}
    yRNAs = {}
    for file in samples:
        f = file.split('/')[-3]
        out_di[f], tLENGTH, tYRNAs, tot_counts[f] = read_lengths_dict(file)
        lengths.update(tLENGTH)
        yRNAs.update(tYRNAs)
    return out_di, lengths, yRNAs, tot_counts


def distributions_output(outPath, lengths, yRNAs, out_di, yRNA_counts):
    '''
    pass
    '''
    yRNA_name = {'yRNA1' : 'RNY1',
                 'yRNA2' : 'RNY2',
                 'yRNA3' : 'RNY3',
                 'yRNA4' : 'RNY4',
                 'yRNA5' : 'RNY5'}
    low = int(min(lengths))
    high = int(max(lengths))
    output = '{}/yRNA_distributions.csv'.format(outPath)
    print 'Output saved as: {}'.format(output)
    with open(output, 'w') as f:
        f.write('sample,yRNA,length,len_per,prevalence\n')
        for samp in out_di:
            for y in yRNAs:
                for l in range(low, high + 1):
                    try:
                        len_per = out_di[samp][y][int(l)]
                    except KeyError:
                        len_per = 0
                    f.write('{},{},{},{},{}\n'.format(samp, 
                                                      yRNA_name[y], 
                                                      l, 
                                                      len_per, 
                                                      yRNA_counts[samp][y]))
    return output


def make_R_figures(r_dir, fi):
    '''
    Create a length distribution image from the data.
    Create a yRNA distribution image from the data.
    '''
    os.system('Rscript --vanilla {}/yrnaDistGraph.R {}'.format(r_dir, fi))


def main(outPath, samples):
    out_di, lengths, yRNAs, tot_counts = process_yRNA_files(samples)
    output = distributions_output(outPath, lengths, yRNAs, out_di, tot_counts)
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
