"""
Reference: https://codereview.stackexchange.com/questions/263416/upgma-tree-building-in-python
UPGMA - Unweighted pair group method using arithmetic averages
Takes as input the distance matrix of species as a numpy array
Returns tree either as a dendrogram or in Newick format
"""

from dataclasses import dataclass
from typing import Self, Generator

import numpy as np


@dataclass
class Leaf:
    """Data structure to store leaf of a UPGMA tree"""

    name: str
    up_height: float = 0.0

    def leaves(self) -> Generator[Self, None, None]:
        yield self

    def __len__(self):
        return 1

    def __str__(self):
        return f'{self.name}:{self.up_height}'


@dataclass
class Node:
    """
    Data structure to store OTU of a UPGMA tree
    """

    left: Self | Leaf
    right: Self | Leaf
    up_height: float = 0.0
    down_height: float = 0.0

    @classmethod
    def from_height(cls, left: Self | Leaf, right: Self | Leaf, height: float):
        left.up_height = height if isinstance(left, Leaf) else height - left.down_height
        right.up_height = height if isinstance(right, Leaf) else height - right.down_height
        return cls(left, right, down_height=height)

    def leaves(self) -> Generator[Leaf, None, None]:
        """
        Method to find the taxa under any given node, effectively equivalent to
        finding leaves of a binary tree. Only lists original taxa and not OTUs.
        """
        yield from self.left.leaves()
        yield from self.right.leaves()

    def __len__(self) -> int:
        """
        Method to define len() of a node.

        Returns the number of original taxa under any given node.
        """
        return len(self.left) + len(self.right)

    def __str__(self) -> str:
        """
        Method to give readable print output
        """
        return f"({self.left},{self.right}):{self.up_height}"


def upgma_tree(dist_matrix: np.ndarray, taxa: list[str]) -> Leaf | Node:
    if not taxa:
        raise ValueError('need at least one leaf (taxon)')

    size = len(taxa)
    work_matrix = np.array(dist_matrix, dtype=float)
    np.fill_diagonal(work_matrix, np.inf)

    if work_matrix.shape != (size, size):
        raise ValueError('distance matrix should be squared in the number of taxa')

    # creating node for each taxa
    nodes = list(map(Leaf, taxa))

    while len(nodes) > 1:
        # finding (row, col) of least dist
        least_id = np.unravel_index(work_matrix.argmin(), work_matrix.shape, "C")
        least_dist = work_matrix[least_id]

        # nodes corresponding to (row, col)
        node1, node2 = (nodes[i] for i in least_id)

        # add OTU with node1 and node2 as children. set heights of nodes
        new_node = Node.from_height(node2, node1, least_dist / 2)
        nodes.remove(node1)
        nodes.remove(node2)
        nodes.append(new_node)

        # create new working distance matrix:
        # 1) adding mean of distances as distance to the new node at the bottom row
        row1, row2 = work_matrix[least_id, :]
        len1 = len(node1)
        len2 = len(node2)
        distances_update = ((row1 * len1) + (row2 * len2)) / (len1 + len2)
        work_matrix = np.vstack((work_matrix, distances_update.reshape(1, -1)))
        # 2) adding mean of distances as distance to the new node at the last column
        distances_update = np.vstack((distances_update.reshape(-1, 1), [np.inf]))
        work_matrix = np.hstack((work_matrix, distances_update))
        # 3) removing rows and columns of the selected nodes
        work_matrix = np.delete(np.delete(work_matrix, least_id, axis=0), least_id, axis=1)

    return nodes[-1]


if __name__ == "__main__":
    # data from table 3 of Fitch and Margoliash, Construction of Phylogenetic trees
    taxa = ["Turtle", "Human", "Tuna", "Chicken", "Moth", "Monkey", "Dog"]
    distances = np.array([
        [0, 19, 27, 8, 33, 18, 13],
        [19, 0, 31, 18, 36, 1, 13],
        [27, 31, 0, 26, 41, 32, 29],
        [8, 18, 26, 0, 31, 17, 14],
        [33, 36, 41, 31, 0, 35, 28],
        [18, 1, 32, 17, 35, 0, 12],
        [13, 13, 29, 14, 28, 12, 0],
    ])
    print(upgma_tree(distances, taxa))