#usr/bin/env python

from numpy import *
from TreeClass import TreeClass
import numpy
from ClusteringUtils import *

numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)

def calculate_Q_matrix(matrice, maxVal):
	#To reformat again with numpy special function. This is ugly
	n= matrice.shape[0] #number of node
	Q_matrix=numpy.zeros(shape=matrice.shape)
	numpy.fill_diagonal(Q_matrix, maxVal)
	i=0; j=0
	while(i<n):
		while(j<n):
			if(i!=j):
				Q_matrix[i,j]= (n-2)*matrice[i,j] - (numpy.sum(matrice[i][numpy.where(matrice[i,:]!=matrice[i,i])])) - (numpy.sum(matrice[:,j][numpy.where(matrice[:,j]!=matrice[j,j])]))
			j+=1
		i+=1
		j=0
	return Q_matrix


def paired_node_distance(matrice, smallest_index):
	i, j = smallest_index
	#i, j are the index of the recently joined node
	n= matrice.shape[0]
	#http://en.wikipedia.org/wiki/Neighbor_joining#equation_2
	#distance from the pair members to the new node second term
	x=(numpy.sum(matrice[i][numpy.where(matrice[i,:]!=matrice[i,i])])) - (numpy.sum(matrice[:,j][numpy.where(matrice[:,j]!=matrice[j,j])]))
	dist_i= 0.5*matrice[i, j] +((0.5/(n-2))*(x))
	dist_j= matrice[i,j]-dist_i
	return dist_i, dist_j


def condense_node_order(matrice, smallest_index, node_order):
	"""
	condenses two nodes in node_order based on smallest_index info
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
	dists = paired_node_distance(matrice,smallest_index)
	nodes = [node1,node2]
	pos = [0,1]

	for ind in pos:
		if nodes[ind].get_children():
			pass#n.add_features(Length= d - n.get_child_at().TipLength)
		else:
			nodes[ind].add_features(Length=dists[ind])
		nodes[ind].add_features(TipLength = dists[ind])
	#combine the two nodes into a new TreeNode object
	new_node = TreeClass()
	new_node.add_child(node1)
	new_node.add_child(node2)
	#replace the object at index1 with the combined node
	node_order[index1] = new_node
	#replace the object at index2 with None
	node_order[index2] = None #distance at i=index2 || j=index2 is large_number
	return node_order


def NJ_cluster(matrice, node_order, large_number, nj_depth=None):
	"""
	Node clustering with NJ
	matrice is a numpy array.
	node_order is a list of PhyloNode objects corresponding to the matrice.
	large_number will be assigned to the matrice during the process and
	should be much larger than any value already in the matrice.
	
	WARNING: Changes matrice in-place.
	WARNING: Expects matrice to already have diagonals assigned to large_number
	before this function is called.
	"""

	# this is for a test, should made it into one function with upgma
	num_entries= len(node_order)
	if not nj_depth or nj_depth>(num_entries-1):
		nj_depth = num_entries-1 #default, do all, same as upgma

	tree=None
	smallest_index=[]
	for i in range(nj_depth):
		smallest_index= find_smallest_index(matrice)
		index1, index2= smallest_index
		#we shouldn't have the smallest index on the diagonal
		#if that happen, reset the diagonal to large number
		if index1 == index2:
			matrice[diag([True]*len(matrice))] = large_number
			smallest_index = find_smallest_index(matrice)
		row_order = condense_node_order(matrice, smallest_index, node_order)
		matrice= condense_matrix(matrice, smallest_index, large_number, method='nj')
		tree = node_order[smallest_index[0]]
	return tree, matrice, smallest_index




def test():
	matrice=numpy.array([[1e305,5,9,9,8],[5,1e305,10,10,9],[9,10,1e305,8,7],[9,10,8,1e305,3],[8,9,7,3,1e305]])
	node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
	t= NJ_cluster(matrice, node_order,1e305)
	print t[0]
    #print t.get_ascii(show_internal=True, compact=False, attributes=['Length', 'TipLength', 'name'])


if __name__ == '__main__':
	test()

	"""Q_matrix=calculate_Q_matrix(matrice,1e305)
	s_index=find_smallest_index(Q_matrix)
	a, b = paired_node_distance(matrice,s_index)
	print a, " and ", b
	matrice=del_row_column(matrice,s_index , 1e305, method='nj')
	print matrice

	Q_matrix=calculate_Q_matrix(matrice,1e305)
	s_index=find_smallest_index(Q_matrix)
	a, b = paired_node_distance(matrice,s_index)
	print a, " and ", b
	matrice=del_row_column(matrice,s_index , 1e305, method='nj')
	print matrice"""
