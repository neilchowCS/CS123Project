import numpy as np

from NeighborJoiningV2 import displayGraph, NeighborJoining2

m= np.array([[0,0.189,0.1100,0.1130,0.2150],[0.1890,0,0.1790,0.1920,0.2110], [0.1100,0.1790,0,0.0940,0.2050], [0.1130,0.1920,0.0940,0,0.2140], [0.2150,0.2110,0.2050,0.2140,0]])

displayGraph(NeighborJoining2(m ,5)[1], ["Gorilla","Orangutan","Human","Chimp","Gibbon"])