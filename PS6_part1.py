#!/usr/bin/env python
# coding: utf-8

import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description="A program to k-merize our data set at different sizes for k and plot the result")
    parser.add_argument("-f", "--file", help="use to specify file name", required=True, type = str)
    parser.add_argument("-k", "--kmer", help="use to specify kmer size", type = int, default = 49)
    return parser.parse_args()


args = get_args()               # calls get_args method from above assigns the arguments to args
GIVEN_FILE = args.file          # assigning file name as string to global varible
GIVEN_KMER = args.kmer          # assigning kmer size as int to global variable

KMER_LENGTH = []                # initialize KMER_LENGTH array as global variable to hold all kmer lengths
KMER_COVERAGE = []              # initialize KMER_COVERAGE array as global variable to hold all kmer converges
CONTIG_LENGTH_ARRAY = []        # initialize CONTIG_LENGTH_ARRAY array as global variable to hold all contig lengths

NL = 0                          # Number of lines in file
NUM_CONTIG = 0
TOTAL_NT_PER_LINE = 0           # calculating total nucleotides by the sum of the length of each sequence line
TOTAL_NT_BY_FORMULA = 0         # calculatin total nucleotides by adjusting the kmer length to represent the physical length  (kcnt= len(sequence) - kmer_size + 1) solving for len(sequence) each time

MAX_CONTIG_LENGTH = 0
MEAN_CONTIG_LENGTH = 0.0
CURRENT_CONTIG_LENGTH = 0
MEAN_COVERAGE_DEPTH = 0.0
N50 = 0

with open(GIVEN_FILE, 'r') as contigFile:
    for line in contigFile:
        line = line.strip()
        if(line[0] == '>'):
            temp = re.search(r'.*_([0-9]{2,5}).*_([0-9]+\.[0-9]+)',line)        # tmep.group(1): kmer_length      temp.group(2): kmer_coverage
            KMER_LENGTH.append(temp.group(1))
            KMER_COVERAGE.append(temp.group(2))

            if(CURRENT_CONTIG_LENGTH > MAX_CONTIG_LENGTH):      # keep a running adjustment for MAX_CONTIG_LENGTH
                MAX_CONTIG_LENGTH = CURRENT_CONTIG_LENGTH

            if(CURRENT_CONTIG_LENGTH != 0):                     # keep record of all of the contig lengths
                CONTIG_LENGTH_ARRAY.append(CURRENT_CONTIG_LENGTH)

            CURRENT_CONTIG_LENGTH = 0                           # reset current contig length
            TOTAL_NT_BY_FORMULA += (int(temp.group(1)) + GIVEN_KMER -1)     # using kcnt formula to calculate length of contig
            MEAN_COVERAGE_DEPTH += float(temp.group(2))             # keep running sum of coverage depth
            NUM_CONTIG += 1                                         # keep track of number of contigs seen
        else:
            CURRENT_CONTIG_LENGTH += len(line)              # for all sequence lines keep a running sum of the lengths until the next header is seen ( works for converted and non-converted fa files)    

            TOTAL_NT_PER_LINE += len(line)                  # double checking my kcnt formula is giving me the right answer
        NL += 1

MEAN_CONTIG_LENGTH = TOTAL_NT_BY_FORMULA/NUM_CONTIG
MEAN_COVERAGE_DEPTH = MEAN_COVERAGE_DEPTH / NUM_CONTIG

print("KMER SIZE:", GIVEN_KMER)
print("TOTAL LENGTH OF GENOME ASSEMBLY:", TOTAL_NT_BY_FORMULA)
print("NUM CONTIG:",NUM_CONTIG)
print("MAX CONTIG LENGTH:",MAX_CONTIG_LENGTH)
print("MEAN CONTIG LENGTH:",MEAN_CONTIG_LENGTH)
print("MEAN COVERAGE DEPTH:",MEAN_COVERAGE_DEPTH)

CONTIG_LENGTH_ARRAY.sort()
temp_n50_sum = 0
for element in CONTIG_LENGTH_ARRAY:             # calculating n50 which is just the weighted median
    temp_n50_sum += element
    if(temp_n50_sum > (TOTAL_NT_BY_FORMULA/2)):
        N50 = element
        break

print("N50: ", N50)

CONTIG_DISTRIBUTION = {}
for element in CONTIG_LENGTH_ARRAY:         # put contig lengths into buckets
    temp_var = element // 100               # doing integer division to simplify the numbers so it is easier to put into buckets
    if(temp_var in CONTIG_DISTRIBUTION.keys()):
        CONTIG_DISTRIBUTION[temp_var] += 1
    else:
        CONTIG_DISTRIBUTION[temp_var] = 1


print('#Contig length\tNumber of contigs in this category')
for index,key in enumerate(sorted(CONTIG_DISTRIBUTION.keys())):
    print(str(key*100) + "\t" + str(CONTIG_DISTRIBUTION[key]))
