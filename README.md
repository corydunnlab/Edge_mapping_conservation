# Edge_mapping_conservation

----
#### Commercial use without the permission of the software authors (Cory Dunn and Bala Ani Akpinar) is forbidden. Academic use is permitted, with citation of the following work:

#### _B.A. Akpınar, Paul O. Carlson, and C.D. Dunn. A novel phylogenetic analysis combined with a machine learning approach predicts human mitochondrial variant pathogenicity. bioRxiv. doi:._
----

**General Information**

**replace-nonstandard.py** – Replaces non-standard nucleotides in a FASTA alignment file with standard nucleotides. The replacement is picked at random, from the IUPAC correspondents that are represented in the alignment for the given position. 

Dependency: Bio.SeqIO

Requirement: “IUPAC_nonstandard_codes.txt” file must be in the same folder where the script is run. The use of different code files is possible. Please refer to source code to get more information.

```
$ python replace-nonstandard.py -a [alignment_file]
```

**extract-correspondence-for-merged-alignment-v.1.1.py** - Generates a correspondence file that reports alignment position correspondences of each aligned residue of a species of interest. This correspondence file is later required as input for "AA-presence-lookup-table-v.1.1.py" and "direct-subst-lookup-table-proteins-v.1.1.py" scripts.

Since this script is designed for merged alignments, (1) a list of protein order in the merged alignment file, and
(2) a FASTA file of reference sequences for the species of interest are required. The FASTA file of reference sequences does not have to be in the same order as given in the list file, but the names should match. If the alignment file contains only one protein, these two files are still required. Example input and output files are provided in the "Examples" folder.

```
$ python extract-correspondence-for-merged-alignment-v.1.1.py -a [alignment_file] -r [reference_FASTA_file] -l [list_file]
-conv <species_convention> -o <output_prefix>
```
Dependency: Bio.SeqIO




 

**AA-presence-lookup-table-v.1.1.py** - Please refer to source code to get more information.

**direct-subst-lookup-table-proteins-v.1.1.py** - Please refer to source code to get more information.

**Usage**

Please refer to source codes to get more information.
Use --help” or “-h” flags to get usage information. 

Contact Bala Ani Akpinar (@aniakpinar) for code-related issues and Cory Dunn (@corydunnlab) for more general inquiries.
