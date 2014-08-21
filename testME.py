from Multipolysolver import *
from TreeLib import *
from TreeLib.SupportUtils import *
import os, shutil

def runTest(output, basedir, align_type, smap, specietree, alignfile, mltree, phylogeny, slimit=-1, plimit=-1, reroot='best', seuil=95, seq=1):

	mltree_ext="%s.root.bootstrap.tree"%align_type

	shutil.copy(os.path.join(basedir, 'RAxML_bipartitions.bootstrap%s.tree'%(align_type)),mltree)
	cmd="raxmlHPC-SSE3 -f I -m GTRGAMMA -t %s -n %s -w %s" %(mltree, mltree_ext[1:], os.path.abspath(basedir))
	executeCMD(cmd)
	shutil.copy(phylogeny, phylogeny.replace(".tree", ".true.tree"))
	phylogeny=phylogeny.replace(".tree", ".true.tree")
	rooted_tree=os.path.join(basedir, os.path.basename(phylogeny).replace(".true.tree", mltree_ext))
	shutil.copy(os.path.join(basedir,"RAxML_rootedTree%s.root.bootstrap.tree"%align_type), rooted_tree)
	
	#mltree="exp/0.align.root.bootstrap.tree"
	gtree=os.path.join(basedir, os.path.basename(phylogeny).replace(".tree", "%s.%s.bootstrap.tree"%(align_type,seuil)))
	maxLTree= TreeClass(mltree)
	maxLTree.contract_tree(seuil=seuil)
	maxLTree.write(outfile=gtree, format=0)

	methodCompare(output, rooted_tree, smap, specietree, alignfile, gtree, seuil, mltree_ext, reroot, slimit, plimit, phylogeny)
	#treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel --show-help
	#treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel --show-help


if __name__ == '__main__':
	
	import argparse

	reroot_option= ['none', 'all', 'best'];
	parser = argparse.ArgumentParser(description='Simulation script, run polytomySolver,TreeFix and compare the output. The input file name should respect TreeFix standard')
	parser.add_argument('-s', '--species', dest='specietree', help="Name of the file containing the species newick tree.",required=True)
	parser.add_argument('-S', '--smap', dest='smap', help="Gene to species map. Use the standard format.",required=True)
	parser.add_argument('-g', '--mltree', dest='mltree', help="Name of the file containing the gene newick tree.")
	parser.add_argument('-t', '--phylotree', dest='realtree', help="Name of the file containing the true phylogenetic tree.")
	parser.add_argument('-a', '--align', dest='alignment', help="Name of the multi-alignment file")
	parser.add_argument('-r', '--reroot', choices=reroot_option, dest='reroot', default='none', help= "Root settings for PolytomySolver")
	parser.add_argument('--slimit', type=int, dest="sol_limit", default=30, help="slimit setting for polytomysolver")
	parser.add_argument('--plimit', type=int, default=-1,  dest="path_limit", help="plimit setting for polytomysolver")
	parser.add_argument('-w', '--workdir', dest='workdir', default="./", help="Working directory. The directory which contains each simulation")
	parser.add_argument('--seq', dest='sequence', type=int, nargs=2, help="In case of multiple directory in the working directory, choose the directory to process.")
	parser.add_argument('-n', dest='align_type', help="Alignement file extension", required=True)
	parser.add_argument('--seuil', dest='seuil', type=int, default=95, help="Support contraction threshold for polytomySolver")

	args = parser.parse_args()

	output=os.path.join(args.workdir, "output.txt")
	def_tree="RAxML_bipartitions.bootstrap%s.tree"%(args.align_type)

	if(args.sequence):
		for seq in xrange(args.sequence[0], args.sequence[1]):
			basedir=os.path.join(args.workdir, str(seq))
			mltree= os.path.join(basedir, "%s%s.bootstrap.tree"%(seq,args.align_type))
			alignfile=os.path.join(basedir, "%s%s"%(seq,args.align_type))
			phylogeny=os.path.join(basedir, "%s.tree"%(seq))
			if(os.path.exists(os.path.join(basedir, def_tree))):
				runTest(output, basedir, args.align_type, args.smap, args.specietree, alignfile, mltree, phylogeny, reroot=args.reroot, seuil=args.seuil, plimit=args.path_limit, slimit=args.sol_limit, seq=seq)

	else:
		runTest(output,args.workdir,args.align_type, args.smap, args.specietree, os.path.join(args.workdir,args.alignment), os.path.join(args.workdir,args.mltree), os.path.join(args.workdir,args.realtree), reroot=args.reroot, seuil=args.seuil, plimit=args.path_limit, slimit=args.sol_limit)
