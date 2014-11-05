#usr/bin/env python

from numpy import *
from TreeClass import TreeClass
import TreeUtils
import numpy
import random as rnd
from Multpolysolver import *
from ClusteringUtils import *

numpy.set_printoptions(precision=3)
numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)

def test():

	names=['c_1', 'd_1', 'd_2','c_2', 'c_3', 'd_3' ,'a_1', 'b_1','a_2', 'b_2', 'e_1', 'e_2', 'f_1', 'f_2', 'g_1', 'g_2', 'h_1', 'h_2', 'i_1', 'i_2']

	for i in range(5):
		node_order = rnd.sample(names, rnd.randint(15,20))
		matrix= makeFakeDstMatrice(len(node_order),0,1,0)
		genetree= TreeUtils.make_random_tree(node_order, contract_seuil=0.5)
		sp_tree= make_random_specie_tree(genetree, sep='_', capitalize=False)
		saveMatrix("testinput/matrix_%i.dist"%i, matrix,node_order)
		genetree.write(outfile="testinput/genetree_%i.nwk"%i)
		sp_tree.write(outfile="testinput/specietree%i.nwk"%i)
		print "\n***------------------------------------------------------------------***"
		print i, ")"
		print "***------------------------------------------------------------------***\n"
		print genetree, "\nMatrix :", matrix, "\nSpecie tree", sp_tree
		nj_sol = solvePolytomy("testinput/genetree_%i.nwk"%i, "testinput/specietree%i.nwk"%i, "testinput/matrix_%i.dist"%i, poly_limit=-1, method='nj', genepos="prefix")
		upgma_sol = solvePolytomy("testinput/genetree_%i.nwk"%i, "testinput/specietree%i.nwk"%i, "testinput/matrix_%i.dist"%i, poly_limit=-1, method='upgma', genepos="prefix")
		t=0
		for sol in nj_sol:
			sol.write(outfile="testinput/nj_sol%i_%i.nw"%(i,t))
			t+=1
		
		t=0
		for sol in upgma_sol:
			sol.write(outfile="testinput/upgma_sol%i_%i.nw"%(i,t))
			t+=1

def make_random_specie_tree(genetree, sep='_', capitalize=False):
	specie_list= list(set(map(lambda x: x.split(sep)[0].upper() if capitalize else x.split(sep)[0] ,genetree.get_leaf_names())))
	sp_tree= TreeClass()
	sp_tree.populate(len(specie_list), specie_list)
	return sp_tree

if __name__ == '__main__':
	test()