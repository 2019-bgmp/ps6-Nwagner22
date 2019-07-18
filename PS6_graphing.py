#!/usr/bin/env python
# coding: utf-8

import argparse
import matplotlib.pyplot as plt
import re

def get_args():
    parser = argparse.ArgumentParser(description="A program to k-merize our data set at different sizes for k and plot the result")
    parser.add_argument("-f", "--file", help="use to specify file name", required=True, type = str)
    return parser.parse_args()


args = get_args()               # calls get_args method from above assigns the arguments to args
GIVEN_FILE = args.file          # assigning file name as string to global varible

BUCKET_SIZES = []
NUM_CONTIGS_IN_BUCKET = []
BUCKET_TITLE = '#Contig length'
NUM_CONTIGS_TITLE = 'Number of contigs in this Bucket'

with open(GIVEN_FILE, 'r') as slurmFile:
    while True:
        line = slurmFile.readline().strip()
        values = line.split()
        if(line == ""):
            break
        if(values[0] == '#Contig'):
            break

    # while True:                                 # run this while loop if the slurm output file has multiple outputs in the same file. Skips the data I have already graphed
    #     line = slurmFile.readline().strip()
    #     values = line.split()
    #     if(line == ""):
    #         break
    #     if(values[0] == '#Contig'):
    #         break

    # while True:                                 # run this while loop if the slurm output file has multiple outputs in the same file. Skips the data I have already graphed
    #     line = slurmFile.readline().strip()
    #     values = line.split()
    #     if(line == ""):
    #         break
    #     if(values[0] == '#Contig'):
    #         break

    while True:
        line = slurmFile.readline().strip()
        values = line.split()
        if(line == ""):
            break
        if(len(values) != 2):
            break
        else:
            BUCKET_SIZES.append(int(values[0]))
            NUM_CONTIGS_IN_BUCKET.append(int(values[1]))


plt.plot(BUCKET_SIZES,NUM_CONTIGS_IN_BUCKET)
plt.xlabel(BUCKET_TITLE)
plt.ylabel(NUM_CONTIGS_TITLE)
plt.title('ins_length:500 cov_cutoff:auto: Distribution of contig lengths')
plt.yscale('log')

plt.savefig('ins_500_cov_auto_distribution.jpg')