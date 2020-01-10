#!/usr/bin/python
# coding=utf-8
#
# B.A. Akpınar†, Paul O. Carlson, and C.D. Dunn†. A novel phylogenetic analysis combined with a machine learning
# approach predicts human mitochondrial variant pathogenicity. bioRxiv. doi:

# This script generates a look-up table for the presence/absence of the reference amino acid among a list of characters
# from descendant nodes and the Root for a each position, in a given alignment.

# Reference amino acid refers to the amino acid of a species of interest at the given position.

# The script requires the direct output of "extract-correspondence-for-merged-alignment-v.1.1.py" script,
# which includes the reference amino acid information (-corr parameter).

# The descendant characters -among which the reference amino acid is searched- may or may not include extant species,
# depending on the -m parameter (internal/all). Default value is "all", meaning that the extant species
# will be included.
# The Root is always included.

# Dependencies: Bio.SeqIO and Bio.Phylo, xlsxwriter

from Bio import SeqIO
from Bio import Phylo
import argparse
import sys
import xlsxwriter

parser = argparse.ArgumentParser(prog='Amino acid presence lookup table',
                                 usage='$python AA-presence-lookup-table-v1.1.py  -a [alignment_file] -t [tree_file] '
                                       '-corr [correspondence_file] -m <internal/all>')

parser.add_argument('--alignment', '-a',
                    help='Please specify an alignment file in FASTA format')
parser.add_argument('--tree', '-t',
                    help='Please specify a tree file in NEWICK format')
parser.add_argument('--correspondence', '-corr',
                    help='Please specify a correspondence file for a species of interest')
parser.add_argument('--mode', '-m', default='all',
                    help='Please specify mode of the search in tree edges: internal/all (default=all)')

args = parser.parse_args()

alignment_file = args.alignment
tree_file = args.tree
correspondence_file = args.correspondence
edge_mode = args.mode


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


# To define ALL edges in a given tree
def all_edges(tree):

    alledges = []
    for parent_node in tree.find_clades(terminal=False, order='level'):
        for child_node in parent_node.clades:
            alledges.append(parent_node.name + '*' + child_node.name)

    return alledges


# To define all INTERNAL edges in a given tree
def all_internal_edges(tree):

    allINTedges = []
    for parent_node in tree.find_clades(terminal=False, order='level'):
        for child_node in parent_node.clades:
            if '#' in child_node.name:
                allINTedges.append(parent_node.name + '*' + child_node.name)

    return allINTedges


# A list of all amino acids

AA_list_str = 'FWMLIVACSTYGPQNHKRDEX'
AA_list = list(AA_list_str)

##########################################################################################

all_values = read_fasta(alignment_file)
alignment_length = len(all_values[all_values.keys()[0]])

# Define the root
my_tree = Phylo.read(tree_file, 'newick')
root = my_tree.common_ancestor(my_tree.find_clades(order='level'))

if edge_mode == 'all':  # If ALL edges are requested
    edges_all = all_edges(my_tree)
    print 'ALL ~descendant~ edges (and the Root) will be included in the computation.'
elif edge_mode == 'internal':  # If only INTERNAL edges are requested
    edges_all = all_internal_edges(my_tree)
    print 'Only INTERNAL ~descendant~ edges (and the Root) will be included in the computation.'
else:
    print 'Unidentified mode (-m) of search for edges.'
    print 'Please check your script parameters.'
    sys.exit()

# Collect all descendant and Root characters for each alignment position
characters_dict = {}
for pos in range(0, alignment_length):
    pos_name = 'Pos' + str(pos + 1)
    characters_dict[pos_name] = [all_values[root.name][pos]]
    for edge in edges_all:
        descendant = edge.split('*')[1]
        characters_dict[pos_name].append(all_values[descendant][pos])

data_all = {}
order_list = []
print_list = []

# Compute presence/absence
with open(correspondence_file) as cfile:
    for line in cfile:
        elems = line.strip('\n').strip('\r').split('\t')
        identifier = elems[2].replace('Homo_sapiens_', '')  # COI_M1
        identifier_short = identifier.split('_')[0]  # COI
        # Printing progress ...
        if identifier_short not in print_list:
            print 'Processing: Homo_sapiens_%s' % identifier_short
            print_list.append(identifier_short)

        alignment_pos = elems[0]  # Pos1
        char_list = characters_dict[alignment_pos]
        AA_info = [identifier]

        for AA in AA_list:
            if AA in char_list:
                AA_info.append(1)  # Present
            else:
                AA_info.append(0)  # Absent

        data_all[identifier] = AA_info
        order_list.append(identifier)

# Writing the output
print '...'
print 'Writing results...'

if edge_mode == 'all':  # If ALL edges are requested
    workbook = xlsxwriter.Workbook('%s_AApresence_table_for_ALLedges.xlsx' % alignment_file.split('.')[0])
else:  # If only INTERNAL edges are requested
    workbook = xlsxwriter.Workbook('%s_AApresence_table_for_INTERNALedges.xlsx' % alignment_file.split('.')[0])

cell_format = workbook.add_format({'font': 'Helvetica Neue'})
row_names = ['\t'] + AA_list

output = []
for item in order_list:
    gene_name = item.split('_')[0]
    if gene_name not in output:
        worksheet = workbook.add_worksheet('%s' % gene_name)
        i = 1
        worksheet.write_column(0, 0, row_names, cell_format)
        worksheet.write_column(0, i, data_all[item], cell_format)
        output.append(gene_name)
        i += 1
    else:
        worksheet.write_column(0, i, data_all[item], cell_format)
        i += 1

workbook.close()
