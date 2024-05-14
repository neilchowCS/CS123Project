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

# Example HLA sequences 
allele_to_seq = {
    'A*01:01': 'ATCGTACGATCG',
    'A*02:01': 'AGTCTGACAGTC',
    'A*03:01': 'ACGTACGTAGCT',
    'A*24:02': 'ACTGACTGACGT',
    'B*08:01': 'GTCAGTCAGTCA',
    'B*07:02': 'GACTGACTGACT',
    'B*15:01': 'GCATGCATGCAT',
    'B*08:03': 'GCTAGCTAGCTA',
    'C*07:01': 'TAGCTAGCTAGC',
    'C*07:02': 'TACGTACGTACG',
    'C*03:03': 'TCAGTCAGTCAG',
    'C*07:04': 'TGACTGACTGAC',
}

# Process each pair of individuals
pairs = hla_data['id'].unique()
for pair_id in pairs:
    pair_data = hla_data[hla_data['id'] == pair_id]
    seq_records = []
    for index, row in pair_data.iterrows():
        seq = (allele_to_seq.get(row['hla_a'], '') +
               allele_to_seq.get(row['hla_b'], '') +
               allele_to_seq.get(row['hla_c'], ''))
        if seq:
            seq_records.append(SeqRecord(Seq(seq), id=row['person']))

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
    print(f"\nPerformance for Pair ID {pair_id}")
    print(performance_df)

    # Print the results of UPGMA and NJ algorithms for analysis
    print(f"\nUPGMA Result for Pair ID {pair_id}:")
    print(upgma_result)

    print(f"\nNeighbor Joining (NJ) Result for Pair ID {pair_id}:")
    print(nj_result)

    # Function to check compatibility based on distance
    def check_compatibility(distance_matrix, threshold=0.1):
        num_seqs = len(distance_matrix)
        compatible_pairs = []
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                if distance_matrix[i, j] <= threshold:
                    compatible_pairs.append((i, j))
        return compatible_pairs

    # Define a compatibility threshold 
    compatibility_threshold = 0.2  

    # Check compatibility based on the distance matrix
    compatible_pairs = check_compatibility(distance_matrix, compatibility_threshold)
    print(f"\nCompatible Pairs for Pair ID {pair_id} (based on threshold):")
    print(compatible_pairs)

    # Visualization helper function for trees
    def visualize_tree(result, title):
        fig, ax = plt.subplots(figsize=(10, 5))
        for i, cluster in enumerate(result):
            ax.plot([i] * len(cluster), cluster, 'o')
        ax.set_title(title)
        plt.show()

    # Visualize UPGMA and NJ results
    visualize_tree(upgma_result, f'UPGMA Clusters for Pair ID {pair_id}')
    visualize_tree(nj_result, f'NJ Clusters for Pair ID {pair_id}')
