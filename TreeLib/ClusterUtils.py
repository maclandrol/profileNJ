#usr/bin/env python

from numpy import *
from TreeClass import TreeClass
import numpy

numpy.set_printoptions(precision=3)
numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)


def find_smallest_index(matrice):
	"""Return smallest number i,j index in a matrice
	A Tuple (i,j) is returned. 
	Warning, the diagonal should have the largest number so it will never be choose
	"""
	#winner = numpy.argwhere(matrice==numpy.amin(matrice)).tolist()
	#for m in winner:
	#	m_reverse=m[:]
	#	m_reverse.reverse()
	#	if(m_reverse in winner):
	#		winner.remove(m)
	
	#print "number of possible solution here :" , len(winner) , "solution :", winner
	return numpy.unravel_index(matrice.argmin(), matrice.shape)


def condense_matrix(matrice, smallest_index, large_value, method='upgma'):
	"""Matrice condensation in the next iteration

	Smallest index is returned from find_smallest_index.
	For both leaf at i and j a new distance is calculated from the average of the corresponding
	row or the corresponding columns
	We then replace the first index (row and column) by the average vector obtained
	and the second index by an array with large numbers so that
	it is never chosen again with find_smallest_index.
	Now the new regroupement distance value is at the first position! (on row and column)
	"""
	first_index, second_index = smallest_index
	#get the rows and make a new vector by updating distance
	rows = take(matrice, smallest_index, 0)

	#default we use upgma
	if(method.lower()=='nj'):
		new_vector = (numpy.sum(rows, 0)-matrice[first_index,second_index])*0.5

	else:
		new_vector = average(rows, 0)
	
	#replace info in the row and column for first index with new_vector
	matrice[first_index] = new_vector
	matrice[:, first_index] = new_vector
	#replace the info in the row and column for the second index with 
	#high numbers so that it is ignored
	matrice[second_index] = large_value
	matrice[:, second_index] = large_value
	return matrice


def del_row_column(matrice, last_fus_index, large_value, method='upgma'):
	"""This function condense a matrix, actualise the distance between node and remove
	node already merged row/column from the matrix"""
	first_index, second_index = last_fus_index
	matrice = condense_matrix(matrice, last_fus_index, large_value, method)
	x = matrice[:-1,:-1]
	x[:second_index,second_index:] = matrice[:second_index,second_index+1:]
	x[second_index:,:second_index] = matrice[second_index+1:,:second_index]
	x[second_index:,second_index:] = matrice[second_index+1:,second_index+1:]
	return x

def remove_ij(x, i, j):
	# Row i and column j divide the array into 4 quadrants
	y = x[:-1,:-1]
	y[:i,j:] = x[:i,j+1:]
	y[i:,:j] = x[i+1:,:j]
	y[i:,j:] = x[i+1:,j+1:]
	return y

def calculate_Q_matrix(matrice, maxVal):
	
	n= matrice.shape[0]
	Q_matrix=numpy.zeros(shape=matrice.shape)
	numpy.fill_diagonal(Q_matrix, maxVal)
	diag_value= numpy.unique(matrice[diag([True]*len(matrice))])
	infinite_pos=numpy.where(numpy.apply_along_axis(lambda x : (x==maxVal).all(), 1, matrice))[0]
	m= infinite_pos.tolist()
	i=0; j=0
	while(i<n):
		while(j<n):
			if(i!=j and not (i in m or j in m) ) :# and not ((matrice[i]==maxVal).all() or (matrice[:,j]==maxVal).all())):
				Q_matrix[i,j]= (n-len(m)-2)*matrice[i,j] - (numpy.sum(matrice[i][numpy.invert(numpy.in1d(matrice[i], diag_value))])) - (numpy.sum(matrice[:,j][numpy.invert(numpy.in1d(matrice[:,j], diag_value))]))
			j+=1
		i+=1
		j=0
	Q_matrix[infinite_pos]=maxVal
	Q_matrix[:,infinite_pos]=maxVal

	return Q_matrix


def paired_node_distance(matrice, smallest_index,maxVal):
	i, j = smallest_index
	#i, j are the index of the recently joined node
	n= matrice.shape[0]
	infinite_pos=numpy.where(numpy.apply_along_axis(lambda x : (x==maxVal).all(), 1, matrice))[0]
	m= infinite_pos.tolist()
	diag_value= numpy.unique(matrice[diag([True]*len(matrice))])
	#http://en.wikipedia.org/wiki/Neighbor_joining#equation_2
	#distance from the pair members to the new node second term
	x=numpy.sum(matrice[i][numpy.invert(numpy.in1d(matrice[i], diag_value))]) - numpy.sum(matrice[:,j][numpy.invert(numpy.in1d(matrice[:,j], diag_value))])
	if(n-len(m)-2>0):
		dist_i= 0.5*matrice[i, j] +((0.5/(n-len(m)-2))*(x))
		dist_j= matrice[i,j]-dist_i
		return dist_i, dist_j
	else:
		#We have only two node to join (final join)
		# Find the index of the node not already joined
		possible_val= range(0, n)
		pair_dist_ind = list(set(possible_val)- set(m)) #the index of the final pair to join
		distance = matrice[pair_dist_ind[0], pair_dist_ind[1]]
		#In this case, we split the dist value by two 
		return distance/2.0, distance/2.0


def condense_node_order(matrice, smallest_index, node_order, method='upgma'):
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
	
	if(method.lower()=='nj'):
		dist = paired_node_distance(matrice,smallest_index,1e305)
	
	else:
		distance = matrice[index1, index2]
		dist= (distance/2.0 , distance/2.0)
	
	nodes = [node1,node2]
	pos = [0,1]

	for ind in pos:
		if nodes[ind].is_leaf():
			nodes[ind].add_features(Length=dist[ind])
		nodes[ind].add_features(TipLength = dist[ind])
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
		Q_matrix=calculate_Q_matrix(matrice,large_number)
		smallest_index= find_smallest_index(Q_matrix)
		index1, index2= smallest_index
		#we shouldn't have the smallest index on the diagonal
		#if that happen, reset the diagonal to large number
		if index1 == index2:
			matrice[diag([True]*len(matrice))] = large_number
			smallest_index = find_smallest_index(matrice)
		row_order = condense_node_order(matrice, smallest_index, node_order, method='nj')
		matrice= condense_matrix(matrice, smallest_index, large_number, method='nj')
		tree = node_order[smallest_index[0]]
	return tree, matrice, smallest_index


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
		row_order = condense_node_order(matrice, smallest_index, node_order, method='upgma')
		matrice = condense_matrix(matrice, smallest_index, large_number, method='upgma')
		tree = node_order[smallest_index[0]]

	return tree, matrice, smallest_index


def treeCluster(matrice, node_order, large_number, depth=None, method='upgma'):
	if(method.lower()=='nj'):
		return NJ_cluster(matrice, node_order, large_number, nj_depth=depth)
	else :
		return UPGMA_cluster(matrice, node_order, large_number, upgma_depth=depth)


def distMatProcessor(dist_file, maxValue, nFlag=False):
	"""Formating distance matrix from a file input and node order for 
		UPGMA join
	"""

	read_fl=False
	dist_matrix= []
	node_order=[]
	with open(dist_file, 'r') as infile:
		x_ind=0
		for line in infile:
			line= line.strip()
			if(line):
				if not read_fl:
					read_fl=True
				else:
					x_ind+=1
					line_list= [getIntValue(x.strip(), x_ind, y_ind, maxValue, nFlag) for y_ind, x in enumerate(line.split())]
					dist_matrix.append(line_list[1:])
					node_order.append(line_list[0])
	
	return numpy.array(dist_matrix), node_order


def makeFakeDstMatrice(n, dmin, dmax, maxVal):
	"""Create a fake distance matrice"""
	b = (dmax-dmin)*numpy.random.random_sample(size=(n,n))+dmin
	b_sym=(b + b.T)/2
	numpy.fill_diagonal(b_sym, maxVal)
	return b_sym


def saveMatrix(filename, matrix, node_order):
	#matrix[numpy.where(matrix==1e305)]=0
	with open(filename, 'w+') as out:
		out.write("\t%i\n"%len(node_order))
		lines=[]
		for entry in matrix.tolist():
			line=node_order.pop(0)+"\t"+ " ".join(map(str,entry))+"\n"
			lines.append(line)
		out.writelines(lines)
	return True


def getIntValue(number, x_ind, y_ind, maxValue, nFlag=False):
	"""Get a distance matrice validate input from a string"""
	try:
		n=float(number)
		return maxValue if ((n<0 and nFlag) or (x_ind==y_ind)) else n
	except ValueError:
		return number
