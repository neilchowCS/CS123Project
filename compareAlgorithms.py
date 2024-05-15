from Bio import SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from matplotlib import pyplot as plt

import NeighborJoiningV2
import UPGMA
from GenerateMatrix import calculate_distance_matrix

sequences = []
names = []

with open("input.txt") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(record)
        names.append(record.id)

alignments = MultipleSeqAlignment(sequences)
distance_matrix = calculate_distance_matrix(alignments)

upgma_result = UPGMA.UPGMA(distance_matrix,names)
nj_result = NeighborJoiningV2.NeighborJoining2(distance_matrix,5)

parser = Phylo.NewickIO.Parser.from_string(upgma_result)
tree = parser.parse()
Phylo.draw(list(tree)[0], branch_labels=lambda c: c.branch_length, do_show=False)
plt.savefig('upgma.png')  # Adjust dpi as needed

NeighborJoiningV2.displayGraph(nj_result[1], names)