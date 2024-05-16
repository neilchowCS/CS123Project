import numpy as np
from Bio.Align import MultipleSeqAlignment

def calculate_distance_matrix(alignments):
    """Generate distance matrix and return as 2D numpy array."""
    num_seqs = len(alignments)
    distance_matrix = np.zeros((num_seqs, num_seqs))
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = str(alignments[i].seq)
            seq2 = str(alignments[j].seq)
            differences = sum(1 for a, b in zip(seq1, seq2) if a != b)
            distance_matrix[i, j] = distance_matrix[j, i] = differences / len(seq1)
    return distance_matrix
