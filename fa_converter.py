#!/usr/bin/env python
# coding: utf-8

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program to convert FASTA files into two lines per record")
    parser.add_argument("-i", "--input", help="use to specify input file name", required=True, type = str)
    parser.add_argument("-o", "--output", help="use to specify output file name", required=True, type = str)    
    return parser.parse_args()


args = get_args()  
INPUT_FILE = args.input          # assigning input file name as string to global varible
OUTPUT_FILE = args.output        # assigning output file name as string to global varible


with open(OUTPUT_FILE, 'w') as outputFile:
    with open(INPUT_FILE, 'r') as inputFile:
        firstLine = inputFile.readline().strip()
        outputFile.write(firstLine + "\n")
        for line in inputFile:
            line = line.strip()
            if(line[0] == '>'):
                outputFile.write("\n" + line + "\n")
            else:
                outputFile.write(line)

