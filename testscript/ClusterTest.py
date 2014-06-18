#usr/bin/env python

from numpy import *
from TreeClass import TreeClass
import numpy
from ClusteringUtils import *

numerictypes = numpy.core.numerictypes.sctype2char
Float = numerictypes(float)


def test():
	matrice=numpy.array([[1e305,5,9,9,8],[5,1e305,10,10,9],[9,10,1e305,8,7],[9,10,8,1e305,3],[8,9,7,3,1e305]])
    node_order=[TreeClass("a:1;"),TreeClass("b:1;"),TreeClass("c:1;"),TreeClass("d:1;"),TreeClass("e:1;") ]
    t_upgma= UPGMA_cluster(matrice, node_order,1e305)
   	t_nj= NJ_cluster(matrice, node_order,1e305)
   	print "***** UPGMA *****\n"
    print t_upgma[0].get_ascii(show_internal=True, compact=False, attributes=['Length', 'TipLength', 'name'])
    print "\n***** NJ *****\n"
    print t_nj[0].get_ascii(show_internal=True, compact=False, attributes=['Length', 'TipLength', 'name'])
    #print t.get_ascii(show_internal=True, compact=False, attributes=['Length', 'TipLength', 'name'])

if __name__ == '__main__':
	test()