from PolytomySolver.Multipolysolver import *
from TreeLib import *
from TreeLib.SupportUtils import *
import os, shutil, glob
from itertools import chain

def runTest(outfile, basedir, align_type, smap, specietree, alignfile, mltree, phylogeny, slimit=-1, plimit=-1, reroot='best', seuil=[100], algo='nj'):
	cleanup(basedir)
	mltree_ext="%s.root.bootstrap.tree"%align_type
	rootreefile="RAxML_rootedTree%s.root.bootstrap.tree"%align_type
	rootreefileinfo="RAxML_info%s.root.bootstrap.tree"%align_type
	try:
		os.remove(os.path.join(basedir,rootreefile))

	except OSError:
		pass

	try:
		os.remove(os.path.join(basedir,rootreefileinfo))

	except OSError:
		pass

	shutil.copy(os.path.join(basedir, 'RAxML_bipartitions.bootstrap%s.tree'%(align_type)),mltree)
	cmd="raxmlHPC-SSE3 -f I -m GTRGAMMA -t %s -n %s -w %s" %(mltree, mltree_ext[1:], os.path.abspath(basedir))
	executeCMD(cmd)
	#shutil.copy(phylogeny, phylogeny.replace(".tree", ".true.tree"))
	#phylogeny=phylogeny.replace(".tree", ".true.tree")
	rooted_tree=os.path.join(basedir, os.path.basename(phylogeny).replace(".tree", mltree_ext))
	shutil.move(os.path.join(basedir,rootreefile), rooted_tree)
	
	#mltree="exp/0.align.bootstrap.tree"
	#rooted_tree = "exp/0.align.root.bootstrap.tree"
	alpha = 95
	maxLikeTree= TreeClass(mltree)
	data_size=len(maxLikeTree.get_leaves())

	with open("%s_treefix.csv"%outfile, 'a+') as OUT:
		header=['tree', 'leaves', 'TruePhylo_dlc','Raxml_dlc', 'Raxml_rf', 'Raxml_maxrf', 'Treefix_dlc', 'Treefix_rf', 'Treefix_maxrf', 'Treefix_dinl', 'Treefix_p-val', 'Treefix_time', 'TruePhylo_nad', 'TruePhylo_dup', 'TruePhylo_lost','Raxml_nad', 'Raxml_dup', 'Raxml_lost','Treefix_nad', 'Treefix_dup', 'Treefix_lost']
		headerWrite(OUT, header)
		treefixtree=treeCompute(OUT, rooted_tree, smap, specietree, alignfile, mltree_ext, phylogeny, data_size, mltree)

	header=['tree', 'PolySolver_dlc','PolySolver_nad','PolySolver_dup','PolySolver_lost', 'PolySolver_rf','PolySolver_maxrf', 'PolySolver_dinl','PolySolver_p-val',  'PolySolver_time', 'ML_time', 'total_polytomy', 'polysolver_nsol', 'best_sol_aLRT', 'is_bestrf',  'best_rf']

	for thres in seuil:
		gtree=os.path.join(basedir, os.path.basename(phylogeny).replace(".tree", "%s.%s.bootstrap.tree"%(align_type,thres)))
		if(thres>100):
			maxLikeTree.toStar()
		else:
			maxLikeTree.contract_tree(seuil=thres)
		polytomy_number=len(list(maxLikeTree.iter_polytomies()))
		maxLikeTree.write(outfile=gtree, format=0)

		with open("%s%s_polytomysolver.csv"%(outfile, thres), 'a+') as PSOUT:
			headerWrite(PSOUT, header)
			pnjbesttree, basename=pscompute(PSOUT, rooted_tree, smap, specietree, alignfile, gtree, thres, mltree_ext, reroot, slimit, plimit, phylogeny, polytomy_number, algo)
		
		if(thres==alpha):
			lkheader=['tree', 'Raxml_logL', 'TreeFix_logL', 'ProfileNJ%s_logL'%(alpha), 'TruePhylo_logL','Raxml_AU', 'TreeFix_AU', 'ProfileNJ%s_AU'%(alpha), 'TruePhylo_AU','Raxml_SH', 'TreeFix_SH', 'ProfileNJ%s_SH'%(alpha), 'TruePhylo_SH']
			with open("%s_likelihood.csv"%outfile, 'a+') as LIKEOUT:
				headerWrite(LIKEOUT, lkheader)
				treesfile=[rooted_tree, phylogeny.replace('.tree', '.treefix.tree') , pnjbesttree, phylogeny]
				time, consel_out=  runRaxmlPval(basedir, alignfile, 4, out=None, listfile=treesfile, sort=1)
				to_write=[basename]
				to_write.extend(list(chain.from_iterable([consel_out['likelihood'], consel_out['au'], consel_out['sh']])))
				LIKEOUT.write("\t".join([str(val) for val in to_write]) + "\n")

	#treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel --show-help
	#treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel --show-help
	#raxmlHPC-SSE3 -f I -m GTRGAMMA -t 0.bootstrap.align.tree -n broot

		
def headerWrite(stream, header):
	line ="\t".join(header)+"\n"
	if line not in stream.readline():
		stream.write(line)


def cleanup(basedir):
	deletethis = []
	print "Cleaning process ..."
	deletethis.extend(glob.glob("%s/*polytomysolver*"%basedir))
	deletethis.extend(glob.glob("%s/*output*"%basedir))
	deletethis.extend(glob.glob("%s/"%basedir +('[0-9]*')+"*.align.*.bootstrap.tree"))
	deletethis.extend(glob.glob("%s/*[pt][stf].ml*"%basedir))
	deletethis.extend(glob.glob("%s/*phy*"%basedir))
	deletethis.extend(glob.glob("%s/*treefix*"%basedir))
	deletethis.extend(glob.glob("%s/*rmt*"%basedir))
	deletethis.extend(glob.glob("%s/*pv*"%basedir))
	deletethis.extend(glob.glob("%s/*perSiteLLs*"%basedir))
	deletethis.extend(glob.glob("%s/*trees"%basedir))

	for file in set(deletethis):
		os.remove(file)
	


def pscompute(outfile, mltree, smap, specietree, alignfile, gtree, seuil, mltree_ext, r_option, slimit, plimit, correctPhylo, polytomy_number, algo):
	w_dir=os.path.dirname(mltree)
	basename, align_ext=name_extractor(os.path.basename(alignfile), ext=".")
	distmat= getDistMatrix(alignfile,os.path.join(w_dir, basename+".dist"))
	ps_out=os.path.join(w_dir, basename+".%s.polytomysolver.tree"%seuil)
	ps_time, rst=runPolytomySolver(gtree, smap, specietree, ps_out, distmat, r_option, slimit, plimit, algo)
	
	treesfile= os.path.join(w_dir, "polysolver_ml.tree")
	n_psol, m_cost=fix_ps_out(ps_out,treesfile)
	#lik_time, likelihoods = phymlikelihood(alignfile, align_ext, ps_out, n_psol)  #!! edit this to go back to phyml version
	
	###Add choice with consel
	phy_seq=os.path.join(w_dir, basename+".phy")
	
	#consel_output=runConsel(w_dir, phy_seq,treesfile, n_psol) #!! edit this for phyml version
	lik_time, consel_output = runRaxmlPval(w_dir, alignfile, n_psol, out=treesfile, sort=9)
	likelihoods=consel_output['likelihood']
	bestposition=0
	psrfs=[]
	psdl_costs=[]
	psnad_costs=[]
	psdup_costs=[]
	pslost_costs=[]
	psdinls=[]
	pspvals=[]
	psmax_rfs=[]
	# compute pval and Dlnl for true tree using RAxML tree to optimize

	with open(os.path.join(w_dir, basename+".%s.polytomysolver.all"%seuil), 'w') as ps_outfile:
		header=['#','dl-cost', 'nad', 'ad', 'lost', 'likelihood', 'rf', 'max_rf', 'dinl', 'pval', 'consel_rank', 'consel_au', 'consel_sh']
		ps_outfile.write("\t".join(header)+"\n")

		for n in xrange(n_psol):
			ps_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s -o %s %s" %("."+align_ext,mltree_ext,".ps.ml%s"%(n+1),".%s.polytomysolver.tree%s"%(seuil, n+1), "%s%s"%(ps_out, (n+1)))
			executeCMD(ps_cmp_lk)
			ps_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".%s.polytomysolver.tree%s"%(seuil,n+1), ".%s.polytomysolver.output%s"%(seuil,n+1), "%s%s"%(ps_out, (n+1)))
			executeCMD(ps_cmp_dl)

			polysolver_rf, polysolver_maxrf=getRFval(correctPhylo, "%s%s"%(ps_out, (n+1)))
			psrfs.append(polysolver_rf); psmax_rfs.append(polysolver_maxrf)
			with open(os.path.join(w_dir, basename+".ps.ml%s"%(n+1)), 'r') as psml, open(os.path.join(w_dir, basename+".%s.polytomysolver.output%s"%(seuil, (n+1))), 'r') as psr:
				pval, dinl=psml.readline().strip().split()
				rdl=psr.readline().strip()
				#assert int(rdl)==m_cost
				psdl_costs.append(rdl); psdinls.append(dinl); pspvals.append(pval)
				ps_nad, ps_dup, ps_lost= retrieveDupAndLostCost("%s%s"%(ps_out, n+1), specietree, smap, sep='_', pos='prefix')
				psnad_costs.append(ps_nad); psdup_costs.append(ps_dup); pslost_costs.append(ps_lost)
				item= int(consel_output['item'][n])-1 if n_psol>1 else 0
				towrite=[(n+1),rdl,ps_nad, ps_dup, ps_lost, likelihoods[n], polysolver_rf, polysolver_maxrf, dinl, pval, consel_output['rank'][item],consel_output['au'][item], consel_output['sh'][item]]
				ps_outfile.write("\t".join([str(val) for val in towrite])+"\n")

	
		#best polysolver result
		bestposition= int(consel_output['item'][0])-1 if n_psol>1 else selectBestTree(likelihoods)
		best_rf=min(psrfs)
		polysolver_dlc=psdl_costs[bestposition]
		polysolver_rf=psrfs[bestposition]
		polysolver_maxrf=psmax_rfs[bestposition]
		polysolver_dinl=psdinls[bestposition]
		polysolver_pval=pspvals[bestposition]
		line= [basename, polysolver_dlc, psnad_costs[bestposition],psdup_costs[bestposition],pslost_costs[bestposition], polysolver_rf, polysolver_maxrf, polysolver_dinl,polysolver_pval, ps_time, lik_time, polytomy_number, n_psol, likelihoods[bestposition], polysolver_rf==best_rf, best_rf]
		outfile.write("\t".join([str(val) for val in line])+"\n")

	return "%s%s"%(ps_out, bestposition+1), basename

def treeCompute(outfile, mltree, smap, specietree, alignfile, mltree_ext, correctPhylo, datasize, unroot_ml_tree):
	
	w_dir=os.path.dirname(mltree)
	basename, align_ext=name_extractor(os.path.basename(alignfile), ext=".")
	logfile= os.path.join(w_dir, basename+".treefix.log")
	tf_out=os.path.join(w_dir, basename+".treefix.tree")
	tf_time, rst=runTreeFix(mltree, smap, specietree, "."+align_ext, mltree_ext, logfile)
	
	#compute likelihood
	tf_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s -o %s %s" %("."+align_ext, mltree_ext, ".tf.ml", ".treefix.tree", tf_out)
	executeCMD(tf_cmp_lk)

	# compute dl cost
	raxml_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, mltree_ext, ".raxml.output", mltree)
	executeCMD(raxml_cmp_dl)

	tf_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".treefix.tree", ".treefix.output", tf_out)
	executeCMD(tf_cmp_dl)

	tt_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".tree", ".true.output", correctPhylo)
	executeCMD(tt_cmp_dl)

	#Custom dup, lost, nad
	tf_nad, tf_dup, tf_lost= retrieveDupAndLostCost(tf_out, specietree, smap, sep='_', pos='prefix')
	tt_nad, tt_dup, tt_lost= retrieveDupAndLostCost(correctPhylo, specietree, smap, sep='_', pos='prefix')
	raxml_nad, raxml_dup, raxml_lost= retrieveDupAndLostCost(mltree, specietree, smap, sep='_', pos='prefix')

	
	with open(os.path.join(w_dir, basename+".tf.ml"), 'r') as tfml,  open(os.path.join(w_dir, basename+".treefix.output"), 'r') as tfr, open(os.path.join(w_dir, basename+".true.output"), 'r') as ttr, open(os.path.join(w_dir, basename+".raxml.output"), 'r') as rxr:
		#ml score
		treefix_pval, treefix_dinl=tfml.readline().strip().split()
		#rf score
		raxml_rf, raxml_maxrf= getRFval(correctPhylo, unroot_ml_tree, unroot=True)
		treefix_rf, treefix_maxrf=getRFval(correctPhylo, tf_out)
		#reconciliation output
		treefix_rdl=tfr.readline().strip()
		raxml_rdl=rxr.readline().strip()
		default_rdl=ttr.readline().strip()
		line=[basename, datasize, default_rdl, raxml_rdl, raxml_rf, raxml_maxrf, treefix_rdl, treefix_rf, treefix_maxrf,treefix_dinl, treefix_pval, tf_time, tt_nad, tt_dup, tt_lost,raxml_nad, raxml_dup, raxml_lost,tf_nad, tf_dup, tf_lost]
		outfile.write("\t".join([str(val) for val in line])+"\n")


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
	parser.add_argument('--seqlist', dest='seqlist', nargs=2, help="Choose the file listing your input and the number of input to process.")
	parser.add_argument('--seuil', dest='seuil', type=int, nargs='+', help="Support contraction threshold for polytomySolver")    
	parser.add_argument('--out', dest='out', default="output.csv", help="Output file of the analysis.")
	parser.add_argument('-c', '--cluster', dest='algo', default="nj", help="Cluster algorithm")
	args = parser.parse_args()

	outfile=os.path.join(args.workdir, args.out)
	def_tree="RAxML_bipartitions.bootstrap%s.tree"%(args.align_type)
	seq_list=[]
	if(args.seqlist):
		instream=open(args.seqlist[0], 'r').read().strip().split()
		seq_list=[int(x) for x in instream][:int(args.seqlist[1])]
	
	elif(args.sequence):
		seq_list=xrange(args.sequence[0], args.sequence[1])

	if(seq_list):
		for seq in seq_list:
			basedir=os.path.join(args.workdir, str(seq))
			mltree= os.path.join(basedir, "%s%s.bootstrap.tree"%(seq,args.align_type))
			alignfile=os.path.join(basedir, "%s%s"%(seq,args.align_type))
			phylogeny=os.path.join(basedir, "%s.tree"%(seq))
			if(os.path.exists(os.path.join(basedir, def_tree))):
				runTest(outfile, basedir, args.align_type, args.smap, args.specietree, alignfile, mltree, phylogeny, reroot=args.reroot, seuil=args.seuil, plimit=args.path_limit, slimit=args.sol_limit, algo=args.algo)
	else:
		runTest(outfile, args.workdir, args.align_type, args.smap, args.specietree, os.path.join(args.workdir,args.alignment), os.path.join(args.workdir,args.mltree), os.path.join(args.workdir,args.realtree), reroot=args.reroot, seuil=args.seuil, plimit=args.path_limit, slimit=args.sol_limit, algo=args.algo)
