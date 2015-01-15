from PolytomySolver import *
from TreeLib import TreeUtils, TreeClass, params
import operator

oritree, specietree, distance_matrix, node_order = TreeUtils.polySolverPreprocessing("tests/testinput/nostar_genetree.tree", "tests/testinput/speciestree.newick", None, capitalize=True)
tree_list=[oritree]
tree_list.extend(oritree.reroot())
dl_costs=[]
for genetree in tree_list:
    dl_cost=Multipolysolver.computePolytomyReconCost(genetree, specietree, verbose=False)
    dl_costs.append(dl_cost)
    print dl_cost

best_dl = min(dl_costs)
print "And our winner is %s, count=%s"%(best_dl, len([x for x in xrange(len(dl_costs)) if dl_costs[x]==best_dl]))
genetree= TreeClass("tests/testinput/result.nw", format=0)
specietree= TreeClass("tests/testinput/speciestree.newick")
