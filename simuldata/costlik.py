#!/usr/bin/env python
import argparse
import re
from TreeLib import TreeClass, TreeUtils

parser = argparse.ArgumentParser(description='cost-likelihood')
parser.add_argument('-t', '--trees', type=argparse.FileType('r'),  dest='trees',
                    help="a file containing all the tree.", required=True)

parser.add_argument('-g', '--sMap', type=argparse.FileType('r'),
                    dest='smap', help="Gene to species map. Use the standard format.", required=True)

parser.add_argument('-s', '--species', dest='species', help="Species tree in newick.", required=True)

parser.add_argument('-l', '-c', '--lik', type=argparse.FileType('r'),
                    dest='consel', help="consel likelihood file", required=True)

args = parser.parse_args()

genetrees = []
ad = []
nad = []
loss = []
specietree = TreeClass(args.species)
speciemap = {}

if(args.smap):
	regexmap = {}
	with (open(args.smap, 'rU') if isinstance(args.smap, basestring) else args.smap) as INPUT:
		for line in INPUT:
			g, s = line.strip().split()
			if ('*') in g and '.*' not in g:
				g = g.replace('*', '.*')
			g_regex = re.compile(g, re.IGNORECASE)
			regexmap[g_regex] = s

for line in args.trees:
	line = line.strip()
	if line.startswith('>'):
		pass
	else :
		genetrees.append(line)
		t = TreeClass(line)
		for leaf in t:
			for key, value in regexmap.iteritems():
				if key.match(leaf.name):
					speciemap[leaf.name] = value
		t.set_species(speciesMap=speciemap, sep=None, speclist=specietree.get_leaf_name())

		mapping = TreeUtils.lcaMapping(t, specietree)
		TreeUtils.reconcile(t, mapping, lost=True)
		n, a, l = TreeUtils.detComputeDupLostScore(t)
		ad.append(a)
		nad.append(n)
		loss.append(l)

consel_out = []
ind = 0
for line in args.consel:
	line = line.strip()
	if ind < 2:
		pass
	elif 'rank' in line:
		line +="\tAppD\tNad\tLosses\tCost\n"
		consel_out.append(line)

	elif line:
		l_list =  line.split()
		items = int(l_list[2])
		line += "\t%2d\t%2d\t%2d\t%2d\n"%(ad[items-1], nad[items-1], loss[items-1], ad[items-1]+ nad[items-1]+ loss[items-1])
		consel_out.append(line)

	ind += 1

for line in consel_out:
	print line