from Multipolysolver import *
from TreeLib import *
from TreeLib.SupportUtils import *
import os, shutil

def runTest(outfile, basedir, align_type, smap, specietree, alignfile, mltree, phylogeny, slimit=-1, plimit=-1, reroot='best', seuil=95):

	mltree_ext="%s.root.bootstrap.tree"%align_type

	shutil.copy(os.path.join(basedir, 'RAxML_bipartitions.bootstrap%s.tree'%(align_type)),mltree)
	cmd="raxmlHPC-SSE3 -f I -m GTRGAMMA -t %s -n %s -w %s" %(mltree, mltree_ext[1:], os.path.abspath(basedir))
	executeCMD(cmd)
	#shutil.copy(phylogeny, phylogeny.replace(".tree", ".true.tree"))
	#phylogeny=phylogeny.replace(".tree", ".true.tree")
	rooted_tree=os.path.join(basedir, os.path.basename(phylogeny).replace(".tree", mltree_ext))
	shutil.copy(os.path.join(basedir,"RAxML_rootedTree%s.root.bootstrap.tree"%align_type), rooted_tree)
	
	#mltree="exp/0.align.root.bootstrap.tree"
	gtree=os.path.join(basedir, os.path.basename(phylogeny).replace(".tree", "%s.%s.bootstrap.tree"%(align_type,seuil)))
	maxLTree= TreeClass(mltree)
	data_size=len(maxLTree.get_leaves())
	maxLTree.contract_tree(seuil=seuil)
	polytomy_number=len(list(maxLTree.iter_polytomies()))

	maxLTree.write(outfile=gtree, format=0)
	methodCompare(outfile, rooted_tree, smap, specietree, alignfile, gtree, seuil, mltree_ext, reroot, slimit, plimit, phylogeny, data_size, mltree, polytomy_number)
	#treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel --show-help
	#treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel --show-help
	#raxmlHPC-SSE3 -f I -m GTRGAMMA -t 0.bootstrap.align.tree -n broot


def methodCompare(outfile, mltree, smap, specietree, alignfile, gtree, seuil, mltree_ext, r_option, slimit, plimit, correctPhylo, datasize, unrootmltree, polytomy_number):
    
        w_dir=os.path.dirname(mltree)
        basename, align_ext=name_extractor(os.path.basename(alignfile), ext=".")
        distmat= getDistMatrix(alignfile,os.path.join(w_dir, basename+".dist"))
        logfile= os.path.join(w_dir, basename+".treefix.log")
        ps_out=os.path.join(w_dir, basename+".polytomysolver.tree")
        tf_out=os.path.join(w_dir, basename+".treefix.tree")

        ps_time=runPolytomySolver(gtree, smap, specietree, ps_out, distmat, r_option, slimit, plimit)
        tf_time=runTreeFix(mltree, smap, specietree, "."+align_ext, mltree_ext, logfile)
        n_psol=fix_ps_out(ps_out)
        likelihoods = phymlikelihood(alignfile, align_ext, ps_out, n_psol)
        psrfs=[]
        psdl_costs=[]
        psdinls=[]
        pspvals=[]
        psmax_rfs=[]
        # compute pval and Dlnl for true tree using RAxML tree to optimize
        tf_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s -o %s %s" %("."+align_ext, mltree_ext, ".tf.ml", ".treefix.tree", tf_out)
        executeCMD(tf_cmp_lk)

        tt_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s %s" %("."+align_ext,mltree_ext ,".tt.ml",correctPhylo)
        executeCMD(tt_cmp_lk)

        # compute dl cost
        raxml_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, mltree_ext, ".raxml.output", mltree)
        executeCMD(raxml_cmp_dl)

        tf_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".treefix.tree", ".treefix.output", tf_out)
        executeCMD(tf_cmp_dl)

        tt_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".tree", ".true.output", correctPhylo)
        executeCMD(tt_cmp_dl)


        with open(os.path.join(w_dir, basename+".polytomysolver.all"), 'w') as ps_outfile:
                header=['#','dl-cost', 'dinl', 'pval', 'likelihood', 'rf', 'max_rf']
                ps_outfile.write("\t".join(header)+"\n")
                for n in xrange(n_psol):
                    ps_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s -o %s %s" %("."+align_ext,mltree_ext,".ps.ml%s"%(n+1),".polytomysolver.tree", "%s%s"%(ps_out, (n+1)))
                    executeCMD(ps_cmp_lk)
                    ps_cmp_dl="treefix_compute --type cost -r -m treefix.models.duplossmodel.DupLossModel -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".polytomysolver.tree", ".polytomysolver.output%s"%(n+1), "%s%s"%(ps_out, (n+1)))
                    executeCMD(ps_cmp_dl)
                    polysolver_rf, polysolver_maxrf=getRFval(correctPhylo, "%s%s"%(ps_out, (n+1)))
                    psrfs.append(polysolver_rf); psmax_rfs.append(polysolver_maxrf)

                    with open(os.path.join(w_dir, basename+".ps.m%s"%(n+1)), 'r') as psml, open(os.path.join(w_dir, basename+".polytomysolver.output%s"%(n+1)), 'r') as psr:
                        dinl, pval=psml.readline().strip().split()
                        rdl=psr.readline().strip()
                        psdl_costs.append(rdl); psdinls.append(dinl); pspvals.append(pval)
                        ps_outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%((n+1),rdl, dinl, pval, likelihoods[n], polysolver_rf, polysolver_maxrf))


        with open(os.path.join(w_dir, basename+".tf.ml"), 'r') as tfml,  open(os.path.join(w_dir, basename+".tt.ml"), 'r') as ttml, open(os.path.join(w_dir, basename+".treefix.output"), 'r') as tfr, open(os.path.join(w_dir, basename+".true.output"), 'r') as ttr, open(os.path.join(w_dir, basename+".raxml.output"), 'r') as rxr:
            #ml score
            treefix_dinl, treefix_pval=tfml.readline().strip().split()
            default_dinl, default_pval=ttml.readline().strip().split()
            #rf score
            raxml_rf, raxml_maxrf= getRFval(correctPhylo, mltree, unroot=True)
            treefix_rf, treefix_maxrf=getRFval(correctPhylo, tf_out)
            #reconciliation output
            treefix_rdl=tfr.readline().strip()
            raxml_rdl=rxr.readline().strip()
            default_rdl=ttr.readline().strip()
            #best polysolver result
            bestposition= selectBestTree(likelihoods)
            polysolver_rdl=psdl_costs[bestposition]
            polysolver_rf=psrfs[bestposition]
            polysolver_maxrf=psmax_rfs[bestposition]
            polysolver_dinl=psdinls[bestposition]
            polysolver_pval=pspvals[bestposition]

            line=[basename, datasize, default_rdl, default_dinl, default_pval, raxml_rdl, raxml_rf, raxml_maxrf, treefix_rdl, treefix_rf, treefix_maxrf,treefix_dinl, treefix_pval, tf_time, polysolver_rdl, polysolver_rf, polysolver_maxrf, polysolver_dinl,polysolver_pval, ps_time,polytomy_number, n_psol]
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
	parser.add_argument('--seuil', dest='seuil', type=int, default=95, help="Support contraction threshold for polytomySolver")

	args = parser.parse_args()

	output=os.path.join(args.workdir, "output.txt")
	def_tree="RAxML_bipartitions.bootstrap%s.tree"%(args.align_type)
	with open(output, 'a') as outfile:
		header=['#tree', '#data', 'TruePhylo dlc', 'TruePhylo dinl', 'TruePhylo p-val', 'Raxml dlc', 'Raxml rf', 'Raxml maxrf', 'Treefix dlc', 'Treefix rf', 'Treefix maxrf', 'Treefix dinl', 'Treefix p-val', 'Treefix time', 'PolySolver dlc', 'PolySolver rf','PolySolver maxrf', 'PolySolver dinl', 'PolySolver p-val', 'PolySolver time', '#polytomy', '#polysolver sol']
		outfile.write("\t".join(header)+"\n")

		if(args.sequence):
			for seq in xrange(args.sequence[0], args.sequence[1]):
				basedir=os.path.join(args.workdir, str(seq))
				mltree= os.path.join(basedir, "%s%s.bootstrap.tree"%(seq,args.align_type))
				alignfile=os.path.join(basedir, "%s%s"%(seq,args.align_type))
				phylogeny=os.path.join(basedir, "%s.tree"%(seq))
				if(os.path.exists(os.path.join(basedir, def_tree))):
					runTest(outfile, basedir, args.align_type, args.smap, args.specietree, alignfile, mltree, phylogeny, reroot=args.reroot, seuil=args.seuil, plimit=args.path_limit, slimit=args.sol_limit)
		else:
			runTest(outfile,args.workdir,args.align_type, args.smap, args.specietree, os.path.join(args.workdir,args.alignment), os.path.join(args.workdir,args.mltree), os.path.join(args.workdir,args.realtree), reroot=args.reroot, seuil=args.seuil, plimit=args.path_limit, slimit=args.sol_limit)
