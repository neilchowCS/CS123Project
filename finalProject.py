import time
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
from io import StringIO
import matplotlib.pyplot as plt

# Fake HLA typing data
fake_data = """id,hla_a,hla_b,hla_c
1,A*01:01,B*08:01,C*07:01
2,A*02:01,B*07:02,C*07:02
3,A*03:01,B*15:01,C*03:03
4,A*24:02,B*08:01,C*07:01
5,A*01:01,B*08:03,C*07:04
"""

# Convert the data to a DataFrame
hla_data = pd.read_csv(StringIO(fake_data))

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

# Function to calculate the distance matrix
def calculate_distance_matrix(alignments):
    num_seqs = len(alignments)
    distance_matrix = np.zeros((num_seqs, num_seqs))
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = str(alignments[i].seq)
            seq2 = str(alignments[j].seq)
            differences = sum(1 for a, b in zip(seq1, seq2) if a != b)
            distance_matrix[i, j] = distance_matrix[j, i] = differences / len(seq1)
    return distance_matrix

# UPGMA algorithm implementation
def upgma(distance_matrix):
    clusters = [[i] for i in range(len(distance_matrix))]
    while len(clusters) > 1:
        min_dist = float('inf')
        to_merge = (None, None)
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                dist = sum(distance_matrix[x][y] for x in clusters[i] for y in clusters[j]) / (len(clusters[i]) * len(clusters[j]))
                if dist < min_dist:
                    min_dist = dist
                    to_merge = (i, j)
        new_cluster = clusters[to_merge[0]] + clusters[to_merge[1]]
        clusters = [clusters[i] for i in range(len(clusters)) if i not in to_merge] + [new_cluster]
    return clusters

# NJ algorithm implementation
def nj(distance_matrix):
    n = len(distance_matrix)
    clusters = [[i] for i in range(n)]
    while n > 2:
        q_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                q_matrix[i, j] = q_matrix[j, i] = (n - 2) * distance_matrix[i, j] - sum(distance_matrix[i, :]) - sum(distance_matrix[j, :])
        min_q = np.unravel_index(np.argmin(q_matrix, axis=None), q_matrix.shape)
        i, j = min_q
        new_dist = [(distance_matrix[i, k] + distance_matrix[j, k] - distance_matrix[i, j]) / 2 for k in range(n) if k not in min_q]
        distance_matrix = np.delete(np.delete(distance_matrix, (i, j), axis=0), (i, j), axis=1)
        distance_matrix = np.vstack([distance_matrix, new_dist])
        distance_matrix = np.hstack([distance_matrix, np.append(new_dist, 0).reshape(-1, 1)])
        clusters = [clusters[k] for k in range(len(clusters)) if k not in min_q] + [clusters[i] + clusters[j]]
        n -= 1
    return clusters

# Calculate the distance matrix
distance_matrix = calculate_distance_matrix(alignments)

# Run and time UPGMA
start_time = time.time()
upgma_result = upgma(distance_matrix)
upgma_time = time.time() - start_time

# Run and time NJ
start_time = time.time()
nj_result = nj(distance_matrix)
nj_time = time.time() - start_time

# Performance comparison
performance_data = {
    'Algorithm': ['UPGMA', 'NJ'],
    'Computation Time (s)': [upgma_time, nj_time]
}

performance_df = pd.DataFrame(performance_data)
print(performance_df)

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
import time
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from io import StringIO
import matplotlib.pyplot as plt
from UPGMA import UPGMA
from NeighborJoining import NeighborJoining
from GenerateMatrix import calculate_distance_matrix

# Fake HLA typing data
fake_data = """id,hla_a,hla_b,hla_c
1,A*01:01,B*08:01,C*07:01
2,A*02:01,B*07:02,C*07:02
3,A*03:01,B*15:01,C*03:03
4,A*24:02,B*08:01,C*07:01
5,A*01:01,B*08:03,C*07:04
"""

# Convert the data to a DataFrame
hla_data = pd.read_csv(StringIO(fake_data))

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
upgma_newick = UPGMA(distance_matrix, [str(i) for i in range(len(seq_records))])
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

# Visualization helper function for trees
def visualize_newick(newick, title):
    parser = Phylo.NewickIO.Parser.from_string(newick)
    tree = parser.parse()
    fig, ax = plt.subplots(figsize=(10, 5))
    Phylo.draw(list(tree)[0], branch_labels=lambda c: c.branch_length, do_show=False, axes=ax)
    ax.set_title(title)
    plt.show()

# Visualize UPGMA result
visualize_newick(upgma_newick, 'UPGMA Tree')
