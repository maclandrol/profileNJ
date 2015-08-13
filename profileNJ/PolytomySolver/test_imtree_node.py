#/usr/bin/env python
from TreeLib import TreeClass, TreeUtils, params
from collections import defaultdict as ddict
import numpy as np

class DynPolySolver():
    def __init__(self, genetree, specietree, lcamap):
        self.genetree = genetree
        self.specietree = specietree
        self.lcamap = lcamap
        self.specietree.label_internal_node()
        self.multiplicities = ddict(int)
        #self.reversedmap = TreeUtils.getReverseMap(lcamap, True)
        self._compute_multiplicity()
        # compute image tree        
        self.image_trees = TreeUtils.getImageTreeNode(genetree, specietree, lcamap)
        self.costTable = None

    def _compute_multiplicity(self):
        for g in self.genetree.get_children():
            s = self.lcamap[g]    
            self.multiplicities[s.name] += 1

    def compute_cost_table(self):
        """Compute cost table"""

        child_len = len(self.genetree.get_children())
        self.costTable = [{} for x in xrange(child_len+1)]
        for node in self.image_trees[self.genetree].traverse("postorder"):
            # case of leaf
            dupcost = params.getdup(node)
            losscost = params.getloss(node)
            c_u = 0 
            mult_node = self.multiplicities[node.name]
            if(node.up != None):
                c_u = node.depth - node.up.depth - 1
            for k in xrange(0, child_len+1):
                w = dupcost
                if(mult_node <= k):
                    w = losscost
                w_bar = min(dupcost+losscost, c_u*losscost)
                dt = w_bar*min(k, mult_node) + w*abs(k-mult_node) 

                if node.is_leaf():
                    self.costTable[k][node.name] = dt
            
                elif len(node.get_children())==1:
                    child_list = []
                    for k_p in xrange(0, child_len+1):
                        child_data = self.costTable[k_p].get(node.get_children()[0].name, 0)
                        child_list.append(dt+child_data + k_p*losscost)
                    self.costTable[k][node.name] = min(child_list)

                elif len(node.get_children()) == 2:
                    child_list = []
                    for k_p in xrange(0, child_len+1):
                        child1_data = self.costTable[k_p].get(node.get_children()[0].name, 0)
                        child2_data = self.costTable[k_p].get(node.get_children()[1].name, 0)
                        child_list.append(child1_data + child2_data+dt)
                    self.costTable[k][node.name] = min(child_list)

    def get_cost(self, s, k=0):
        """return cost for a particular specie"""
        try:
            if(s.isinstance(basestring)):
                return self.costTable[k][s]
            else:
                return self.costTable[k][s.name]
        except ValueError:
            print("Can't find cost for %s with %d genes.\n Have you already computed the cost table ?"%(s, k))
            raise
        except:
            raise

    def print_cost_table(self):
        """ print the cost table """
        assert self.costTable != None, "costTable not defined"
        l = len(self.costTable)
        for nodename in self.costTable[0].keys():
            out = nodename + "\t" + "\t".join([str(self.costTable[k][nodename]) for k in xrange(l)])
            print(out)


class LinPolySolver():

    def __init__(self, genetree, specietree, lcamap):
        """initialize just like the super class"""
        self.genetree = genetree
        self.specietree = specietree
        self.lcamap = lcamap
        self.specietree.label_internal_node()
        self.multiplicities = ddict(int)
        #self.reversedmap = TreeUtils.getReverseMap(lcamap, True)
        self._compute_multiplicity()
        # compute image tree        
        self.image_trees = TreeUtils.getImageTreeNode(genetree, specietree, lcamap)
        self.ingene = ddict(int)
        self.outgene = ddict(int)

    def _compute_multiplicity(self):
        for g in self.genetree.get_children():
            s = self.lcamap[g]    
            self.multiplicities[s.name] += 1

    def compute_cost_table(self):
        """compute cost table"""
        alpha = ddict(int)
        beta = ddict(int)
        # first step (1) : traverse Image node in post-order
        for node in self.image_trees[self.genetree].traverse("postorder"):
            mult_node = self.multiplicities[node.name]
            ap, bp, app, bpp = 0, 0, 0, 0
            # node has 2 children
            if(len(node.get_children()) == 2):
                node1 = node.get_children()[0]
                node2 = node.get_children()[1]
                # compute ap and bp for node1
                if(node1.depth - node.depth >= 3):
                    ap = 0; bp = 0;
                elif(node1.depth - node.depth == 2):
                    ap = 0; bp = alpha[node1.name];
                else:
                    ap = alpha[node1.name]; bp = beta[node1.name]

                # compute app and bpp for node2
                if(node2.depth - node.depth >= 3):
                    app = 0; bpp = 0;
                elif(node2.depth - node.depth == 2):
                    app = 0; bpp = alpha[node2.name];
                else:
                    app = alpha[node2.name]; bpp = beta[node2.name]
                alpha[node.name] = mult_node + min(max(ap, app), min(bp, bpp))
                beta[node.name] = mult_node + max(max(ap, app), min(bp, bpp))
            
            # case where node has only one child
            elif(len(node.get_children()) == 1):
                node1 = node.get_children()[0]
                if(node1.depth - node.depth >= 2):
                    alpha[node.name] = mult_node
                    beta[node.name] = mult_node
                else:
                    alpha[node.name] = mult_node
                    beta[node.name] = mult_node + alpha[node1.name]

            # case where node is a leaf 
            else :
                alpha[node.name] = mult_node - 1
                beta[node.name] = mult_node - 1

        # step 2 : traversing in preorder
        # in(u) and out(u) denotes the nos. of genes 
        # flowing into and out of the branch (p(u), u) 
        for node in self.image_trees[self.genetree].traverse("preorder"):
            mult_node = self.multiplicities[node.name]
            if node.is_root():
                self.ingene[node.name] = 0
                self.outgene[node.name] = mult_node + alpha[node.name]
            else :
                p_node = node.up
                if(node.depth - p_node.depth - 1 > 0):
                    self.ingene[node.name] = 0
                else:
                    self.ingene[node.name] = self.outgene[p_node.name] - self.multiplicities[p_node.name]

                self.outgene[node.name] = np.median([self.outgene[p_node.name]-self.multiplicities[p_node.name], alpha[node.name], beta[node.name]])
            self.ingene[node.name] += 1
            self.outgene[node.name] += 1

    def get_in_out(self, s):
        """return cost for a particular specie"""
        try:
            if(s.isinstance(basestring)):
                return self.ingene[s], self.outgene[s]
            else:
                return self.ingene[s.name], self.outgene[s.name]
        except ValueError:
            print("Can't find cost for %s.\n Have you already computed the cost table ?"%(s))
            raise
        except:
            raise

    def print_cost_table(self):
        """ print the cost table """
        assert self.ingene != None and self.outgene!= None, "In and out gene not computed"
        for nodename in self.ingene.keys():
            out = nodename + "\t%d\t%d"%(self.ingene[nodename], self.outgene[nodename])
            print(out)


if __name__ == '__main__':
    
    specnw = "((a,b)e, (c,d)f)g;"
    genenw = "(a_1, a_2, a_3, b_1, (b_2, c_1));"
    genetree = TreeClass(genenw, format=1)
    specietree = TreeClass(specnw, format=1)
    genetree.set_species(sep='_', pos="prefix")

    lcamap = TreeUtils.lcaMapping(genetree, specietree, False)

    dps = DynPolySolver(genetree, specietree, lcamap)
    dps.compute_cost_table()
    dps.print_cost_table()
