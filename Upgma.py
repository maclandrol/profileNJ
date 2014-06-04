#usr/bin/env python

from numpy import *
from TreeClass import TreeClass
import numpy
from ClusteringUtils import *

numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)
   

def condense_node_order(matrice, smallest_index, node_order):
    """condenses two nodes in node_order based on smallest_index info
    
    This function is used to create a tree while condensing a matrice
    with the condense_matrix function. The smallest_index is retrieved
    with find_smallest_index. The first index is replaced with a node object
    that combines the two nodes corresponding to the indices in node order.
    The second index in smallest_index is replaced with None.
    Also sets the branch length of the nodes to 1/2 of the distance between
    the nodes in the matrice"""
    index1, index2 = smallest_index
    node1 = node_order[index1]
    node2 = node_order[index2]
    #get the distance between the nodes and assign 1/2 the distance to the
    #Length property of each node
    distance = matrice[index1, index2]
    nodes = [node1,node2]
    d = distance/2.0
    for n in nodes:
        if n.get_children():
            pass#n.add_features(Length= d - n.get_child_at().TipLength)
        else:
            n.add_features(Length=d)
        n.add_features(TipLength = d)
    #combine the two nodes into a new TreeNode object
    new_node = TreeClass()
    new_node.add_child(node1)
    new_node.add_child(node2)
    #replace the object at index1 with the combined node
    node_order[index1] = new_node
    #replace the object at index2 with None
    node_order[index2] = None #distance at i=index2 || j=index2 is large_number
    return node_order

def UPGMA_cluster(matrice, node_order, large_number, upgma_depth=None):
    """cluster with UPGMA
    
    matrice is a numpy array.
    node_order is a list of PhyloNode objects corresponding to the matrice.
    large_number will be assigned to the matrice during the process and
    should be much larger than any value already in the matrice.
    
    WARNING: Changes matrice in-place.
    WARNING: Expects matrice to already have diagonals assigned to large_number
             before this function is called.
    """
    num_entries = len(node_order)
    if not upgma_depth or upgma_depth>(num_entries-1):
        upgma_depth = num_entries-1 #default, do all
    tree = None
    smallest_index=[]
    for i in range(upgma_depth):
        smallest_index = find_smallest_index(matrice)
        index1, index2 = smallest_index
        #if smallest_index is on the diagonal set the diagonal to large_number
        if index1 == index2:
            matrice[diag([True]*len(matrice))] = large_number
            smallest_index = find_smallest_index(matrice)
        row_order = condense_node_order(matrice, smallest_index, node_order)
        matrice = condense_matrix(matrice, smallest_index, large_number, method='upgma')
        tree = node_order[smallest_index[0]]

    return tree, matrice, smallest_index

def test():
	matrice= numpy.array([[ 1e305, 2, 4, 5, 3 ],[2, 1e305, 4, 4, 4 ], [4., 4, 1e305, 4, 3 ], [ 5, 4, 4, 1e305, 3. ], [ 3, 4, 3, 3, 1e305 ] ])
	node_order=[TreeClass("1:1;"),TreeClass("2:1;"),TreeClass("3:1;"),TreeClass("4:1;"),TreeClass("5:1;") ]
	t= UPGMA_cluster(matrice, node_order,1e305)
	print t
    #print t.get_ascii(show_internal=True, compact=False, attributes=['Length', 'TipLength', 'name'])


if __name__ == '__main__':
	test()