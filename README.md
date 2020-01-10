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
$python replace-nonstandard.py -a [alignment_file]
```

**extract-correspondence-for-merged-alignment-v.1.1.py** - Please refer to source code to get more information.

**AA-presence-lookup-table-v.1.1.py** - Please refer to source code to get more information.

**direct-subst-lookup-table-proteins-v.1.1.py** - Please refer to source code to get more information.

**Usage**

Use --help” or “-h” flags to get usage information. 

Contact Bala Ani Akpinar (@aniakpinar) for code-related issues and Cory Dunn (@corydunnlab) for more general inquiries.
