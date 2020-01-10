#!/usr/bin/python
# coding=utf-8
#
# B.A. Akpınar†, Paul O. Carlson, and C.D. Dunn†. A novel phylogenetic analysis combined with a machine learning
# approach predicts human mitochondrial variant pathogenicity. bioRxiv. doi:

# This script replaces non-standard nucleotides in a FASTA alignment file with standard nucleotides.
# The replacement is picked at random, from the IUPAC correspondents that are represented in the alignment
# for the given position.

# If a position contains only N's and gaps, N's are replaced with gaps.

# Dependency: Bio.SeqIO, IUPAC_nonstandard_codes.txt file

import argparse
from Bio import SeqIO
import random

parser = argparse.ArgumentParser(prog='Replace non-standard',
                                 usage='$python replace-nonstandard.py -a [alignment_file]')

parser.add_argument('--alignment', '-a',
                    help='Please specify an nucleotide alignment file')


args = parser.parse_args()

alignment_file = args.alignment

##########################################################################################
#################### Functions to be used in the script ##################################

# To load the FASTA alignment file


def read_fasta(alignment):  # To read the FASTA alignment file

    aa_dict = {}
    with open(alignment, mode='r') as handle:

        # Using Biopython's parse function to reduce memory footprint
        for record in SeqIO.parse(handle, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            aa_dict[identifier] = sequence

    return aa_dict
##########################################################################################


IUPAC_nonstandard = {}
# If a different code file will be used, the following filename should be edited accordingly.
# The new code files must, however, still be in the same format as the "IUPAC_nonstandard_codes.txt"
# and in the same folder where the script is run.
with open('IUPAC_nonstandard_codes.txt') as codes:
    for line in codes:
        IUPAC_nonstandard[line.split('\t')[0]] = line.strip('\n').strip('\r').split('\t')[1].split(',')
        # IUPAC_nonstandard[Y] = C,T # Example


all_values = read_fasta(alignment_file)
alignment_length = len(all_values[all_values.keys()[0]])

results = open(alignment_file.split('.')[0] + '_StandardOnly.fa', 'w')

for species in all_values.keys():
    reconstructed_sequence = []
    for pos in range(0, alignment_length):  # Process one position at a time
        aa = all_values[species][pos].upper()
        if aa in IUPAC_nonstandard.keys():
            pos_values = []
            for sp in all_values.keys():
                pos_values.append(all_values[sp][pos].upper())

            selection = []  # Make a list of selectable characters found within the tree
            for possible in IUPAC_nonstandard[aa]:
                if possible in pos_values:
                    selection.append(possible)

            if aa == 'N':
                if len(selection) == 0:
                    replacement = '-'  # If the position contains only N's and gaps, replace with gap.
                else:
                    replacement = random.choice(selection)  # Else, randomly pick from the selection list.
            else:
                replacement = random.choice(selection)

            reconstructed_sequence.append(replacement)
        else:
            reconstructed_sequence.append(aa)
    results.write('>' + species + '\n' + ''.join(reconstructed_sequence) + '\n')

results.close()
print 'All done.'
