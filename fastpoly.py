#!/usr/bin/env python
import argparse
import re

from profileNJ.TreeLib import *
from profileNJ.PolytomySolver import *

class GTreeSolver():

    def __init__(self, genetree, specietree, mode="lafond", dupcost=1, losscost=1):

        self.genetree = genetree
        self.specietree = specietree
        self.mode = mode
        self.dupcost = dupcost
        self.losscost = losscost
        if mode == 'zheng':
            specietree.label_internal_node()
            self.lcamap = TreeUtils.lcaMapping(genetree, specietree, False)
            self.solver = SingleSolver.LinPolySolver(genetree, specietree, self.lcamap)
        
        elif mode == 'dynamiq':
            specietree.label_internal_node()
            self.lcamap = TreeUtils.lcaMapping(genetree, specietree, False)
            self.solver = SingleSolver.DynPolySolver(genetree, specietree, self.lcamap, self.dupcost, self.losscost) 
        
        elif mode == "lafond" :
            self.lcamap = TreeUtils.lcaMapping(genetree, specietree, False)
            self.solver = PolySolver.GeneTreeSolver(genetree, specietree, self.lcamap, self.dupcost, self.losscost) 
            self.solver.labelInternalNodes(genetree)
            self.solver.labelInternalNodes(specietree)
            self.solver.use_dp = False
        
        elif mode == "notung":
            specietree.label_internal_node()
            self.lcamap = TreeUtils.lcaMapping(genetree, specietree, False)
            self.solver = SingleSolver.NotungSolver(genetree, specietree, self.lcamap,self.losscost, self.dupcost)

    # to consider everything fair, this is the method
    # I think we should time to get running time for each algorithm
    def solvePolytomies(self, nsol):
        if self.mode == 'lafond':
            return [x+";" for x in self.solver.solvePolytomies(nsol)]
        else:
            # current notung DP implementation does not enable enable 
            # multiple solution.
            return [self.solver.reconstruct()]
        



parser = argparse.ArgumentParser(description='PolySolver : A tool to resolve polytomies in a genetree')
parser.add_argument('-s', '--spectree', dest='specnw',
                    help="Name of the file containing the species newick tree.", required=True)
parser.add_argument('-S', '--sMap', dest='smap', help="Gene to species map. Use the standard format.")
parser.add_argument('-g', '--genetree', dest='genenw', help="Name of the file containing the gene newick tree.", required=True)
parser.add_argument('--sep', dest='gene_sep', help="Specify a gene separator if you're are not using a smap")
parser.add_argument('--losscost', type=float, default=1, dest='losscost', help="Specify the losses cost")
parser.add_argument('--dupcost', type=float, default=1, dest='dupcost', help="Specify the duplication cost")
parser.add_argument('--nsol', type=int,  default=1,  dest='nsol', help="Number of solution to output")
parser.add_argument('--spos', dest='spos', default='prefix', help="Gene position when you have specified a separator. Default value is prefix")
parser.add_argument('-o', '--output', dest='outfile',
                    help="Name of your output files with the corrected tree. The resolutions are printed on stdout if omitted.")
parser.add_argument('--mode', dest='mode', default="lafond", choices=['lafond', 'zheng', 'dynamiq', 'notung'], help="Execution mode of the ")
args = parser.parse_args()

genetree = TreeClass(args.genenw)
specietree = TreeClass(args.specnw)
genetree.set_species()

if(args.smap):
    regexmap = {}
    with open(args.smap, 'rU') as INPUT:
        for line in INPUT:
            g, s = line.strip().split()
            if ('*') in g and '.*' not in g:
                g = g.replace('*', '.*')
            g_regex = re.compile(g, re.IGNORECASE)
            regexmap[g_regex] = s
    for leaf in genetree:
        for key, value in regexmap.iteritems():
            if key.match(leaf.name):
                speciemap[leaf.name] = value

if args.gene_sep:
    genetree.set_species(sep=args.gene_sep, pos=args.spos)

else:
    genetree.set_species(use_fn = lambda x : x.name)

gsolver = GTreeSolver(genetree, specietree, args.mode, args.dupcost, args.losscost)
solutions = gsolver.solvePolytomies(args.nsol)

for sol in solutions:
    tree = TreeClass(sol)
    print tree
    # this is just for testing the complete cost
    #print tree
    tree.set_species(use_fn = lambda x : x.name)
    specietree.lcaprocess = False
    lcamap = TreeUtils.lcaMapping(tree, specietree, multspeciename=False)
    dup, loss = TreeUtils.computeDL(tree, lcamap)
    print "DUPS=",dup
    print "LOSSES=", loss
    print "DL=", dup + loss
