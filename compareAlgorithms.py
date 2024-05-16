import Bio.Align
from Bio import SeqIO, Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from matplotlib import pyplot as plt

import GenerateMatrix
import NeighborJoiningV2
import UPGMA

#suffix of input files read and png files output
suffix = ["HLAA", "HLAB","HLAC"]

for c in range(len(suffix)):
    names = []
    #read sequences in file [c]
    alignment = AlignIO.read("io/aligned_" + suffix[c] + ".txt", "fasta")
    for record in alignment:
        names.append(record.id)

    #calculate distance matrix
    distance_matrix = GenerateMatrix.calculate_distance_matrix(alignment)

    #params: distance matrix, name of each sequence, places to round branch lengths to
    #returns: newick string
    upgma_result = UPGMA.UPGMA(distance_matrix,names, 4)
    # params: distance matrix, places to round branch lengths to
    # returns: tuple (array cluster, networkx graph)
    nj_result = NeighborJoiningV2.NeighborJoining2(distance_matrix,4)

    #save tree to file and queues up plot
    parser = Phylo.NewickIO.Parser.from_string(upgma_result)
    tree = parser.parse()
    Phylo.draw(list(tree)[0], branch_labels=lambda c: c.branch_length, do_show=False)
    plt.savefig('io/upgma_' + suffix[c] + '.png')

    #save tree to file and queues up plot
    NeighborJoiningV2.saveGraph(nj_result[1], names, 'io/nj_' + suffix[c] + '.png')

    #display queued trees
    plt.show()