#usr/bin/python
from numpy import *
from TreeLib import ClusterUtils, TreeClass
import numpy

numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)

def test():
	matrice=numpy.array([[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]])
	node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
	t_upgma= ClusterUtils.UPGMA_cluster(matrice, node_order)
	matrice=numpy.array([[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]])
	node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
	t2_upgma=ClusterUtils.treeCluster(matrice, node_order, depth=None, method='upgma')
	print "***** UPGMA *****\n"
	print t_upgma[0].get_ascii(show_internal=True, compact=False, attributes=['length', 'TipLength', 'name'])
	print t2_upgma[0].get_ascii(show_internal=True, compact=False, attributes=['length', 'TipLength', 'name'])


	matrice=numpy.array([[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]])
	node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
	t_nj = ClusterUtils.NJ_cluster(matrice, node_order)
	matrice=numpy.array([[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]])
	node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
	t2_nj = ClusterUtils.treeCluster(matrice, node_order, depth=None, method='nj')
	print "\n***** NJ *****\n"
	#t_nj[0].unroot()
	print t_nj[0].get_ascii(show_internal=True, compact=False, attributes=['length', 'TipLength', 'name'])
	print t2_nj[0].get_ascii(show_internal=True, compact=False, attributes=['length', 'TipLength', 'name'])


	matrice=numpy.array([[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]])
	node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
	t_rand=ClusterUtils.treeCluster(matrice, node_order, depth=None, method='rand')
	print "\n***** RAND *****\n"
	#t_nj[0].unroot()
	print t_rand[0].get_ascii(show_internal=True, compact=False, attributes=['length', 'TipLength', 'name'])

if __name__ == '__main__':
	test()
