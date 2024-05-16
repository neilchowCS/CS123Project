HLA Sequence Comparison and Analysis

This CS123A project performs a comparison analysis of HLA sequences using traditional alignment methods (UPGMA and Neighbor Joining) to determine which HLA sequences match the most. This analysis can potentially be used to find the most matching organ for people or animals.

Files
compareAlgorithms.py: Main script to perform the analysis.
UPGMA.py: Script containing the implementation of the UPGMA algorithm.
NeighborJoiningV2.py: Script containing the implementation of the Neighbor Joining algorithm.
msa_clustalomega.py: Script containing the Clustal Omega Multiple Sequence Alignment algorithm.
io Folder: contains input and output files
    -input: contains unaligned sequence data in fasta format, input to msa_clustalomega
    -aligned: contains sequence data that is aligned in fasta format, input to compareAlgorithms
    -nj/upgma.png: saved images of phylogenetic trees, output of compareAlgorithms

Prerequisites
Python 3.x
numpy library
pandas library
biopython library
matplotlib library
networkx library
You can install the required Python libraries using pip:
pip install numpy pandas biopython matplotlib networkx

Running the Analysis
1. Enter all HLA A,B,C sequences into their respective input_HLA files.
2. Run msa_clustalomega.py. Aligned sequences will be placed into aligned_HLA files.
If aligned sequences are already available, they can be placed directly into the aligned_HLA files
3. Run compareAlgorithms.py.
This will step through all aligned files, for each one displaying both UPGMA and neighbor joining phylogenetic trees.
Trees will also be stored in their respective .png files.
If you wish to analyze more sequence types or modify which types are analyzes, modify the
suffix = ["HLAA", "HLAB","HLAC"] at the top of compareAlgorithms.
This list represents the suffix of the input files read and tree images returned.
Ex. If you wish to additionally analyze BRAF, you can place it into the list
["HLAA", "HLAB","HLAC", "BRAF"]
and make sure "aligned_BRAF.txt" containing the sequences is in the io folder.

Output
The script will print the performance comparison of UPGMA and Neighbor Joining algorithms and the results for each pair of individuals. Additionally, it will visualize the UPGMA and Neighbor Joining clusters for each pair.