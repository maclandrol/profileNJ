#usr/bin/env python

from numpy import *
from TreeClass import TreeClass
import numpy

numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)


def find_smallest_index(matrice):
	"""Return smallest number i,j index in a matrice
	A Tuple (i,j) is returned. 
	Warning, the diagonal should have the largest number so it will never be choose
	"""
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
	matrice= numpy.delete(matrice, second_index, 0)
	matrice=numpy.delete(matrice,second_index , 1)
	return matrice



def makeFakeDstMatrice(n, dmin, dmax, maxVal):
	"""Create a fake distance matrice"""
	import numpy as np
	b = (dmax-dmin)*np.random.random_sample(size=(n,n))+dmin
	b_sym=(b + b.T)/2
	np.fill_diagonal(b_sym, maxVal)
	return b_sym


def distMatProcessor(dist_file, maxValue):
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
					line_list= [getIntValue(x.strip(), x_ind, y_ind, maxValue) for y_ind, x in enumerate(line.split())]
					dist_matrix.append(line_list[1:])
					node_order.append(line_list[0])
	
	return numpy.array(dist_matrix), node_order


def getIntValue(number, x_ind, y_ind, maxValue):
	"""Get a distance matrice validate input from a string"""
	try:
		n=float(number)
		return n if (n!=0 or x_ind!=y_ind) else maxValue
	except ValueError:
		return number
