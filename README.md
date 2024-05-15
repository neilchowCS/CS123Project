HLA Sequence Comparison and Analysis

This CS123A project performs a comparison analysis of HLA sequences using traditional alignment methods (UPGMA and Neighbor Joining) to determine which HLA sequences match the most. This analysis can potentially be used to find the most matching organ for people or animals.

Files
finalProject.py: Main script to perform the analysis.
UPGMA.py: Script containing the implementation of the UPGMA algorithm.
NeighborJoining.py: Script containing the implementation of the Neighbor Joining algorithm.
GenerateMatrix.py: Script to calculate the distance matrix from sequence alignments.
hla_data.csv: Sample CSV file containing HLA typing data.

Prerequisites
Python 3.x
numpy library
pandas library
biopython library
matplotlib library
You can install the required Python libraries using pip:
pip install numpy pandas biopython matplotlib

Running the Analysis
Prepare HLA Data and HLA Sequences
Run the Main Script: python finalProject.py

Output
The script will print the performance comparison of UPGMA and Neighbor Joining algorithms and the results for each pair of individuals. Additionally, it will visualize the UPGMA and Neighbor Joining clusters for each pair.