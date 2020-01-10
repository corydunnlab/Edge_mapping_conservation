#!/usr/bin/python
# coding=utf-8
#
# B.A. Akpınar†, Paul O. Carlson, and C.D. Dunn†. A novel phylogenetic analysis combined with a machine learning
# approach predicts human mitochondrial variant pathogenicity. bioRxiv. doi:

# This script generates a look-up table for direct substitutions involving the reference amino acid
# for a given position.

# For each position, the script goes through edges in the tree (ALL or only INTERNAL, depending on the -m parameter)
# and searches for edges that include the reference amino acid, in the descendant or ancestor node. When such an edge
# is encountered, the amino acid of the other node is assigned '1' in the look-up table.

# The presence of edges with the reference amino acid on both nodes assigns '1' to the reference amino acid
# in the look-up table.
# Amino acids with '0' in the look-up table indicate that the reference amino acid has never substituted for those
# amino acids in the tree (or vice versa, AA>Reference AA).

# Reference amino acid refers to the amino acid of a species of interest at the given position.

# The script requires the direct output of "extract-correspondence-for-merged-alignment-v.1.1.py" script,
# which includes the reference amino acid information (-corr parameter).


# Dependencies: Bio.SeqIO and Bio.Phylo, xlsxwriter

from Bio import SeqIO
from Bio import Phylo
import argparse
import sys
import xlsxwriter

parser = argparse.ArgumentParser(prog='Direct substitution look-up table for proteins',
                                 usage='$python Direct-subst-lookup-table-proteins-v1.1.py -a [alignment_file] '
                                       '-t [tree_file] -corr [correspondence_file] -m <internal/all>')

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


# To read the FASTA alignment file
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

# Read the FASTA alignment file and record amino acid values for each amino acid position
all_values = read_fasta(alignment_file)

# Read the tree file and retrieve all edges in the tree
my_tree = Phylo.read(tree_file, 'newick')
if edge_mode == 'all':  # If ALL edges are requested
    edges_all = all_edges(my_tree)
    print 'All edges in the tree will be searched for direct substitutions involving the reference amino acid: ' \
          '%r edges' % len(edges_all)
elif edge_mode == 'internal':  # If only INTERNAL edges are requested
    edges_all = all_internal_edges(my_tree)
    print 'Only internal edges in the tree will be searched for direct substitutions involving ' \
          'the reference amino acid: %r edges' % len(edges_all)
else:
    print 'Unidentified mode (-m) of search for edges.'
    print 'Please check your script parameters.'
    sys.exit()

data_all = {}
order_list = []
print_list = []
with open(correspondence_file) as cfile:
    for line in cfile:
        elems = line.strip('\n').strip('\r').split('\t')
        identifier = elems[2].replace('Homo_sapiens_', '')  # COI_M1
        identifier_short = identifier.split('_')[0]  # COI
        ref_aa = list(identifier.split('_')[1])[0]  # M
        alignment_pos = elems[0]  # Pos1
        alignment_pos_numeric = int(alignment_pos.replace('Pos', '')) - 1
        AA_info = [identifier]
        # Printing progress ...
        if identifier_short not in print_list:
            print 'Processing: Homo_sapiens_%s' % identifier_short
            print_list.append(identifier_short)

        for AA in AA_list:
            for edge in edges_all:  # Go through edges in the tree (All/Internal depending on -m)
                parent = edge.split('*')[0]
                child = edge.split('*')[1]
                if all_values[parent][alignment_pos_numeric] == ref_aa:
                    if all_values[child][alignment_pos_numeric] == AA:
                        decision = 1
                        break
                    else:
                        decision = 0
                elif all_values[child][alignment_pos_numeric] == ref_aa:
                    if all_values[parent][alignment_pos_numeric] == AA:
                        decision = 1
                        break
                    else:
                        decision = 0
                else:
                    decision = 0

            AA_info.append(decision)

        data_all[identifier] = AA_info
        order_list.append(identifier)

print '...'
print 'Writing results...'

if edge_mode == 'all':  # If ALL edges are requested
    workbook = xlsxwriter.Workbook('%s_DirectSubst_LookUpTable_ALLedges.xlsx' % alignment_file.split('.')[0])
else:  # If only INTERNAL edges are requested
    workbook = xlsxwriter.Workbook('%s_DirectSubst_LookUpTable_INTedges.xlsx' % alignment_file.split('.')[0])

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
