from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import pandas as pd

# fake HLA typing data
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
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignments)

# Construct the phylogenetic tree using UPGMA algorithm
constructor = DistanceTreeConstructor()
tree = constructor.upgma(distance_matrix)

# Visualize the tree
Phylo.draw(tree)

# Assume traditional HLA typings have a base precision score
traditional_precision_scores = {
    'A': 0.90,
    'B': 0.85,
    'C': 0.80,
}

# Placeholder function for calculating precision improvement from phylogenetic analysis
def calculate_precision_improvement(hla_type):
    # For demonstration purposes, we'll add a flat 0.05 improvement
    # this would be based on the phylogenetic tree analysis later
    return 0.05

# Comparative analysis
for index, row in hla_data.iterrows():
    for locus in ['hla_a', 'hla_b', 'hla_c']:
        # Convert the locus[-1] to uppercase to match the dictionary keys
        traditional_precision = traditional_precision_scores[locus[-1].upper()]
        phylogenetic_improvement = calculate_precision_improvement(row[locus])
        enhanced_precision = traditional_precision + phylogenetic_improvement

        print(f"Sample ID: {row['id']}, Locus: {locus}, Traditional Precision: {traditional_precision}, "
              f"Phylogenetic Precision: {enhanced_precision}")