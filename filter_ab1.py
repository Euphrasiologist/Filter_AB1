#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import os
import sys

# add arguments
parser = argparse.ArgumentParser(description="Process an ab1 or abi file or a directory of files.", 
                                 epilog="Thank you for using!")
parser.add_argument('path', help="Input file(s) path.")
parser.add_argument('-q', '--phred', type=int, help="Phred quality score cutoff, default 5", default=5)
parser.add_argument('-r', '--replacement', type=str, help="The replacement letter for nucleotides below quality threshold. Default=N", default="N")
parser.add_argument('-p', '--proportion', type=float, help="If a sequence has a greater than or equal to this % missing data, sequence will be removed.", default=1)
args = parser.parse_args()


    ######################################
#                 Housekeeping             #
    ######################################

phred_score = args.phred
prop = args.proportion
replace = args.replacement

if prop is not None:
    if prop < 0 or prop > 1:
        sys.exit('Please enter a decimal proportion between zero and one.')

    ######################################
#                 Functions                #
    ######################################

# list to string see https://www.geeksforgeeks.org/python-program-to-convert-a-list-to-string/
def listToString(s):  
    # initialize an empty string 
    str1 = ""  
    # traverse in the string   
    for ele in s:  
        str1 += ele   
    # return string   
    return str1  

# simple fasta printing, may need to be modified if lots of information is wanted in the header.
def print_fasta(seq, id):
    print(">" + id)
    print(seq)

    ######################################
#                 Begin loop                #
    ######################################

# loop over files in the directory
for ab1 in os.listdir(args.path):
    if ab1.endswith(".ab1") or ab1.endswith(".abi"):
        # parse current file
        record = SeqIO.read(os.path.join(args.path, ab1), "abi")
        # which are the bases which are lower than the threshold?
        phred_bool = []
        for phred in record.letter_annotations['phred_quality']:
            phred_bool.append(phred > phred_score)
        # create a new list of nucelotides
        # index will tell us which base we are currently on.
        new_seq = []
        index = 0 
        for nucleotide in record.seq:
            if phred_bool[index] == True:
                new_seq.append(nucleotide)
            else:
                new_seq.append(replace)
            index += 1
        # turn to string
        new_seq2 = listToString(new_seq)
        # calculate proportion missing data
        length = len(new_seq2)
        missing = new_seq2.count(replace)

        if missing/length >= prop:
            continue
        print_fasta(new_seq2, str(record.name + "_" + str(length)))
    
    else:
        sys.exit("No files detected")