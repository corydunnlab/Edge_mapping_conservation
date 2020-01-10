#!/usr/bin/python
# coding=utf-8
#
# B.A. Akpınar†, Paul O. Carlson, and C.D. Dunn†. A novel phylogenetic analysis combined with a machine learning
# approach predicts human mitochondrial variant pathogenicity. bioRxiv. doi:

# This script reports alignment position correspondences of each residue in a protein (or merged proteins)
# of a species of interest.

# Since the script is designed for merged alignments, (1) a list of protein order in the merged alignment file, and
# (2) a FASTA file of reference sequences for the species of interest are required.
# The FASTA file of reference sequences does not have to be in the same order as given in the list file,
# but the names should match.

# If the alignment file contains only one protein, these two files are still required.

# Script usage:
# $ python extract-correspondence-for-merged-alignment-v.1.1.py -a [alignment_file] -r [reference_file] -l [list_file]
# -conv <species_convention> -o <output_prefix>

# Dependency: Bio.SeqIO

import argparse
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(prog='extract correspondence',
                                 usage='$python extract-correspondence-for-merged-alignment-v.1.1.py'
                                       ' -a [alignment_file] -r [reference_file] -l [list_file] '
                                       '-conv <species_convention> -o <output_prefix>')

parser.add_argument('--alignment', '-a',
                    help='Please specify an alignment file in FASTA format')
parser.add_argument('--reference', '-r',
                    help='Please specify a reference file for the merged sequences')
parser.add_argument('--list', '-l',
                    help='Please specify a list file indicating the order of merged sequences')
parser.add_argument('--convention', '-conv',
                    help='Please specify a species of interest for convention')
parser.add_argument('--output', '-o',
                    help='Please specify a short output prefix')

args = parser.parse_args()

alignment_file = args.alignment
ref_file = args.reference
species = args.convention
output_prefix = args.output
lst_file = args.list

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


# Read the FASTA file of reference sequences for the species of interest (convention)
ref_sequences = read_fasta(ref_file)

# Read the merged alignment file
all_values = read_fasta(alignment_file)
species_sequence = all_values[species]
all_values.clear()

ref_dict = {}
p = 1

# Record the order of proteins in the merged alignment file
order_info = []
with open(lst_file) as lst:
    for line in lst:
        order_info.append(line.strip('\n').strip('\r'))

for ref in order_info:
    print 'Reference %s starts at %r' % (ref, p)
    for i, base in enumerate(list(ref_sequences[ref])):
        ref_dict[p] = 'Homo_sapiens_' + ref + '_' + base + str(i + 1)
        p += 1

results = open('%s_%s_correspondence.txt' % (output_prefix, species), 'w')

# Extract correspondences ...

j = 1
for i in range(0, len(species_sequence)):
    if species_sequence[i] != '-':
        align_pos = 'Pos' + str(i + 1)
        species_value = species_sequence[i]
        ref_correspondence = ref_dict[j]
        ref_value = list(ref_correspondence.split('_')[3])[0]
        if species_value == ref_value:
            results.write('\t'.join([align_pos, species_value + str(j), ref_correspondence]) + '\n')
        else:
            print 'Non-matching values with the reference file!', i
            print 'Check your reference sequences.'
            sys.exit()
        j += 1


results.close()
