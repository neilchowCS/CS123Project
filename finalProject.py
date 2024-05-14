import time
import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import matplotlib.pyplot as plt
from UPGMA import UPGMA
from NeighborJoining import NeighborJoining
from GenerateMatrix import calculate_distance_matrix

# Load real HLA typing data from CSV file
hla_data = pd.read_csv('hla_data.csv')

# Map HLA alleles to arbitrary sequences
allele_to_seq = {
    'A*01:01': 'ATCG',
    'A*02:01': 'AGTC',
    'A*03:01': 'ACGT',
    'A*24:02': 'ACTG',
    'B*08:01': 'GTCA',
    'B*07:02': 'GACT',
    'B*15:01': 'GCAT',
    'B*08:03': 'GCTA',
    'C*07:01': 'TAGC',
    'C*07:02': 'TACG',
    'C*03:03': 'TCAG',
    'C*07:04': 'TGAC',
}

# Convert HLA data to SeqRecord objects
seq_records = [SeqRecord(Seq(allele_to_seq[x.hla_a] + allele_to_seq[x.hla_b] + allele_to_seq[x.hla_c]),
                         id=str(index)) for index, x in hla_data.iterrows()]

# Create a MultipleSeqAlignment object
alignments = MultipleSeqAlignment(seq_records)

# Calculate the distance matrix
distance_matrix = calculate_distance_matrix(alignments)

# Run and time UPGMA
start_time = time.time()
upgma_result = UPGMA(distance_matrix, [str(i) for i in range(len(seq_records))])
upgma_time = time.time() - start_time

# Run and time NJ
start_time = time.time()
nj_result = NeighborJoining(distance_matrix)
nj_time = time.time() - start_time

# Performance comparison
performance_data = {
    'Algorithm': ['UPGMA', 'NJ'],
    'Computation Time (s)': [upgma_time, nj_time]
}

performance_df = pd.DataFrame(performance_data)
print(performance_df)

# Print the results of UPGMA and NJ algorithms for analysis
print("\nUPGMA Result:")
print(upgma_result)

print("\nNeighbor Joining (NJ) Result:")
print(nj_result)

# Visualization helper function for trees
def visualize_tree(result, title):
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, cluster in enumerate(result):
        ax.plot([i] * len(cluster), cluster, 'o')
    ax.set_title(title)
    plt.show()

# Visualize UPGMA and NJ results
visualize_tree(upgma_result, 'UPGMA Clusters')
visualize_tree(nj_result, 'NJ Clusters')

# Function to check compatibility based on distance
def check_compatibility(distance_matrix, threshold=0.1):
    num_seqs = len(distance_matrix)
    compatible_pairs = []
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            if distance_matrix[i, j] <= threshold:
                compatible_pairs.append((i, j))
    return compatible_pairs

# Define a compatibility threshold (you can adjust this value)
compatibility_threshold = 0.1

# Check compatibility based on the distance matrix
compatible_pairs = check_compatibility(distance_matrix, compatibility_threshold)
print("\nCompatible Pairs (based on threshold):")
print(compatible_pairs)
