HLA Sequence Comparison and Analysis

This CS123A project performs a comparison analysis of HLA sequences using traditional alignment methods (UPGMA and Neighbor Joining) to determine which HLA sequences match the most. This analysis can potentially be used to find the most matching organ for people or animals.

Files
compareAlgorithms.py: Main script to perform the analysis.
UPGMA.py: Script containing the implementation of the UPGMA algorithm.
NeighborJoiningV2.py: Script containing the implementation of the Neighbor Joining algorithm.
msa_clustalomega.py: Script containing the Clustal Omega Multiple Sequence Alignment algorithm.
GenerateMatrix.py: Script to calculate the distance matrix from sequence alignments.
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

Running Individual Tree Algorithms
GenerateMatrix:
1. Obtain aligned sequences
2. Call the calculate_distance_matrix function
PARAMETERS: Bio.Align.MultipleSeqAlignment
RETURNS: Distance matrix

NeighborJoiningV2:
1. Obtain distance matrix
2. Call the NeighborJoining2 function
PARAMETERS: (matrix: 2D distance matrix, digits_to_round: int)
RETURNS: tuple (clusters: array representation of tree, G: NetworkX graph)
3. Call the saveGraph function to save graph to png
PARAMETERS: (G: NetworkX graph, names: list of taxa names, same count as sequences, file_name: file name to save to)
4. OR Call the displayGraph function to display graph
PARAMETERS: (G: NetworkX graph, names: list of taxa names, same count as sequences)

UPGMA:
1. Obtain distance matrix
2. Call the UPGMA function
PARAMETERS: (matrix: 2D distance matrix, names: list of taxa names, same count as sequences, digits_to_round: int)
RETURNS: string: Newick format representation of tree
3. Use Bio.Phylo to draw tree from Newick
4. Use matplotlib to save graph