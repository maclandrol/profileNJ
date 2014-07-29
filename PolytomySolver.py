#!/usr/bin/env python

import argparse, linecache
from Multipolysolver import *
from TreeLib import *
import sys
from operator import itemgetter

"""
PolytomySolver

Given a non-binary gene tree, a species tree and a gene distance matrix, PolytomySolver outputs all possible
binarization of the gene tree that minimizes the duplication+loss score. The algorithm use clustering algorithm 
like NJ and UPGMA distance criterion to join the 'nearest' subtrees first. If the tree is treated as unrooted (using the -r argument),
the program tries every possible root (this takes a while) in order to find the one that yields the lowest 
DL-score after correction, or to output every rooted correction.

The output file contains the newick tree separated by a new line for all the solutions.

The output format is the following
> Rooted tree 1 ; cost = 
newick solution1
newick solution2

> Rooted tree 2 ; cost = 
newick solution1
newick solution2

"""
class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):

	def _split_lines(self, text, width):
		# this is the RawTextHelpFormatter._split_lines
		if text.startswith('C|'):
			return text[2:].splitlines()  
		return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)

class Output(object):
	def __init__(self, file=None):
		if(file):
			self.out=file
		else:
			self.out= sys.stdout

	def write(self, line):
		self.out.write('%s\n'%line)
	
	def close(self):
		if self.out is not sys.stdout:
			self.out.close()

	@staticmethod
	def error(message):
		sys.stderr.write("Error: %s\n" % message)
		sys.exit(1)

reroot_option= ['none', 'all', 'best'];
parser = argparse.ArgumentParser(description='Polytomy solver with multiple solutions.',formatter_class=SmartFormatter)
parser.add_argument('-s', '--sFile', type=argparse.FileType('r'),  dest='specietree', help="Name of the file containing the species newick tree.",required=True)
parser.add_argument('-g', '--gFile', type=argparse.FileType('r'),  dest='genetree', help="Name of the file containing the gene newick tree.",required=True)
parser.add_argument('-d', '--dist', type=argparse.FileType('r'),  dest='distfile', help="Name of the file containing the distances between each pair of genes (The gene set should be the same for the leaf set of the genetree).",required=True)
parser.add_argument('-o', '--output', type=argparse.FileType('w'), dest='outfile',help="Output file with the corrected tree. The genetree is printed on stdout if omitted.")
parser.add_argument('-gl', '--gLine', type=int,  dest='gline', help="Index of the line in the gene tree file that corresponds to the current gene tree.", default=1)
parser.add_argument('-r', '--reroot', choices=reroot_option, dest='reroot', default='none', help= '''C|Enable/Disable root mode.\n\tnone: disable reroot mode, correct the input polytomies and return the result.\n\tall: enable reroot mode, reroot the genetree at each node and return all polytomy corrected version for each rooted tree.\n\tbest: enable reroot mode, rerrot the genetree at each node and return all polytomy corrected version for the rooted tree with the smallest Dup-Lost score (First one if not unique).\n\n''')
parser.add_argument('-nf', '--nnflag', action='store_true',  dest='nflag', help="Treat negative distances as large distances.")
parser.add_argument('--sep', dest='gene_sep', help="Gene-Specie separator for each leaf name in the genetree. PolytomySolver will guess by default in a very limited list of special character. ***';;' is not a valid separator for the newick format! IF YOUR SEPARATOR IS \";;\", DON'T USE THIS FLAG. THE PROGRAM WILL AUTOMATICALLY GUESS. ***")
parser.add_argument('--mValue', type=float, dest='mval', default=1e305, help="Set largest value in the distance matrix. Entries on the main diagonal and negative values will be replaced by mValue.")
parser.add_argument('-c', '--cluster', choices=['nj', 'upgma'], default='nj', help="C|Set the clustering methods.\n\tupgma: UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering algo.\n\tnj: neighbor joining clustering method, (slower).\n\n")
parser.add_argument('--slimit', type=int, dest="sol_limit", default=30, help="Set the max number of solution per genetree. Possible values are -1 to return all the solution or n, n>0 for a specific number of solution. Setting this argument to -1 is computationally expensive.")
parser.add_argument('--plimit', type=int, default=-1,  dest="path_limit", help="Set the max number of solution for each polytomy in the genetree. Possible values are -1 to explore all the solution or n, n>0 for a specific number of solution. Setting this argument to -1 is also computationally expensive.")
parser.add_argument('-v',  action='store_true', dest='verbose', help=" Output verbosity")
parser.add_argument('--cap', dest='cap',  action='store_true', help="Capitalize the species name of the genetree leaves to match each species. Almost all functions are case sensitive.")

args= parser.parse_args()


genetree = linecache.getline(args.genetree.name, args.gline)
oritree, specietree, distance_matrix, node_order = TreeUtils.polySolverPreprocessing(genetree, args.specietree.name, args.distfile.name, capitalize=args.cap, gene_sep= args.gene_sep, nFlag= args.nflag)
tree_list=[oritree]
outlog= Output(args.outfile)

if args.reroot.lower() == 'all':
	tree_list.extend(oritree.reroot())

elif args.reroot.lower() == 'best':
	tree_list.extend(oritree.reroot())
	dl_costs=[]
	for genetree in tree_list:
		sol= solvePolytomy(genetree, specietree, distance_matrix, node_order, sol_limit=1, method='upgma', path_limit=1, verbose=False)
		dl_costs.append(sol[0].cost)

	best_dl = min(enumerate(dl_costs), key=itemgetter(1))[0]
	tree_list= [tree_list[best_dl]]

count=0
for genetree in tree_list:
	first=True
	count+=1
	polysolution = solvePolytomy(genetree, specietree, distance_matrix, node_order, sol_limit=args.sol_limit, method=args.cluster, path_limit=args.path_limit, verbose= args.verbose, maxVal=args.mval)
	outlog.write('>Tree %s; cost=%s'%(count, polysolution[0].cost))
	for tree in polysolution:
		outlog.write(tree.write())
	
outlog.close()