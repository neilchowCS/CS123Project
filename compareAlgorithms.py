import Bio.Align
from Bio import SeqIO, Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline
from matplotlib import pyplot as plt

import NeighborJoiningV2
import UPGMA
from GenerateMatrix import calculate_distance_matrix

names = []
alignment = AlignIO.read("aligned.txt", "fasta")
for record in alignment:
    names.append(record.id)
    print(record.id)

distance_matrix = calculate_distance_matrix(alignment)

upgma_result = UPGMA.UPGMA(distance_matrix,names, 4)
nj_result = NeighborJoiningV2.NeighborJoining2(distance_matrix,4)

parser = Phylo.NewickIO.Parser.from_string(upgma_result)
tree = parser.parse()
Phylo.draw(list(tree)[0], branch_labels=lambda c: c.branch_length, do_show=False)
plt.savefig('upgma.png')  # Adjust dpi as needed
#plt.close()

NeighborJoiningV2.saveGraph(nj_result[1], names)
plt.show()