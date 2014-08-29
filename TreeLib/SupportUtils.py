import subprocess
import TreeUtils
from TreeClass import TreeClass
import re, sys, os,time
import numpy as np
from Bio.Phylo.Applications import _Phyml
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC
SEQUENCE_ALPHABET = {'dna':IUPAC.unambiguous_dna, 'rna':IUPAC.unambiguous_rna, 'protein':IUPAC.protein}


phymllk='phyml_lk.txt'
phymltree='_phyml_tree.txt'
phymlmodel='_phyml_stats.txt'
phymlbootstrap='_phyml_boot_trees.txt'
phymlbootstrapmodel='_phyml_boot_stats.txt'
phymltrees='_phyml_trees.txt'


def timeit(func):

    def timed(*args, **kw):
        tstart = time.time()
        result = func(*args, **kw)
        tend = time.time()
        ttime= tend-tstart

        print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
        return ttime

    return timed


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

def getDistMatrix(alignfile, outfile):
    cmd="fastdist -I fasta -O phylip -e -o %s %s"%(outfile, alignfile)
    executeCMD(cmd)
    return outfile


def name_extractor(filename, ext="."):
    prefix, ext, suffix = filename.partition(ext)
    return prefix, suffix


def phymlikelihood(sequence, align_ext, treepath, n_sol):
    
    sequence=convert_to_phylip(sequence, 'fasta', 'align')
    likelihoods=[]
    output_stats = "%s_phyml_stats.txt" %sequence
    output_tree = "%s_phyml_tree.txt" %sequence
    ll_keyword = ". Log-likelihood:"
 
    for n in xrange(n_sol):
        phyml = _Phyml.PhymlCommandline(cmd="phyml", input=sequence, input_tree="%s%s"%(treepath,(n+1)) , optimize="none", bootstrap=0)
        phyml()

    with open(output_stats) as search:
        for line in search:
            if ll_keyword in line:
                line = line.replace(ll_keyword, "")
                likelihoods.append(float(line))

    return likelihoods

#raxmlHPC-SSE3 -f I -m GTRGAMMA -t 0.bootstrap.align.tree -n broot

def selectBestTree(likelihoods):
    return likelihoods.index(max(likelihoods))

@timeit
def runPolytomySolver(mltree, smap, spectree, outfile, distmat, r_option, slimit, plimit):
    cmd="python PolytomySolver.py -s %s -S %s -g %s -d %s -o %s -n -r %s --slimit=%s --plimit=%s"%(spectree, smap, mltree, distmat, outfile, r_option, slimit, plimit)
    executeCMD(cmd)


@timeit
def runTreeFix(mltree, smap, spectree, align_ext, mltree_ext, logfile):
    cmd="treefix -s %s -S %s -A '%s' -o %s -n %s -V 2 -l %s %s"%(spectree, smap, align_ext, mltree_ext, ".treefix.tree", logfile, mltree)
    executeCMD(cmd)


def getRFval(refTree_path, tree_path, unroot=False):
    refTree=TreeClass(refTree_path)
    tree=TreeClass(tree_path)
    rf, max_rf, c, p1, p2= refTree.robinson_foulds(tree, unrooted_trees=unroot)

    ###test:
    rooting= tree.reroot()
    for t in rooting:
        rrf, rmax_rf, rc, rp1, rp2= refTree.robinson_foulds(t, unrooted_trees=True)
        if(rrf!=rf):
            print "... rerooting don't yeld the same rf score for tree: %s"%tree_path
            break

    return rf, max_rf


def fix_ps_out(ps_out):
    nsol=0
    with open(ps_out, 'r') as infile:
        for line in infile:
            if(line.startswith('>')):
                pass
            else:
                nsol+=1
                with open("%s%s"%(ps_out, nsol), 'w') as outfile:
                    outfile.write(line)
    return nsol


def write_al_in_nxs_file(fastafile, outnxs="seq.nxs", al=0):

    if al:
        print
        print 'CONVERTING your fasta file to nexus format ... with "clustalw"'
        subprocess.call(['clustalw', ('-INFILE='+fastafile), ('-OUTFILE='+outnxs), '-convert', '-OUTPUT=NEXUS'])

    else:
        print
        print 'Sequence not aligned!!'
        print 'ALIGNING sequences with clustalw ...'
        clustalcmd="clustalw -INFILE=" +fastafile+" -OUTFILE="+outnxs+ " -OUTPUT=NEXUS"
        print clustalcmd
        error= executeCMD(clustalcmd)
        if not error:
            print "DONE! :  sequences ALIGNED with clustalw"
        else:
            print "FAILED at clustalw execution !!!"
            return


    return nexrepair(outnxs)


def convert_to_phylip(filename, filetype, old_ext, rm=False):
    con_file=filename.replace(old_ext, "phy")
    with open(filename, 'rU') as infile, open(con_file, 'w') as outfile:
        al_input= AlignIO.parse(infile,  filetype)
        AlignIO.write(al_input, outfile, 'phylip')
    if(rm):
        os.remove(filename)
    return con_file


def construct_phyml_tree(input_tree):
    phyml= _Phyml.PhymlCommandline(input=input_tree, datatype='nt', bootstrap=100, optimize='tlr')
    print str(phyml)
    phyml()


def executeCMD(cmd):
    p=subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    error=0
    while True:
        out = p.stderr.read(1)
        if p.poll() != None:
            print "Process terminated"
            break

        if out != '':
            sys.stdout.write(out)
            sys.stdout.flush()

        else :
            print "Something went wrong"
            error=1
            break

    return error


def nexus_repair(nexus_file_path):
    with open(nexus_file_path, "r") as nexus_file:
        with open(nexus_file_path+".repaired", "w") as repaired:
            for line in nexus_file:
                #PhyML does not support the "missing" or "gap" format parameters
                if not "missing" or not "gap" in line:
                    repaired.write(line)

    os.remove(nexus_file_path)
    os.rename(nexus_file_path+".repaired", nexus_file_path)
    return nexus_file_path



def clustalo(geneSeq_file_path, treeid, alignment_out_path="", dist_matrix_out_path="", aligned=False, cmd_path="utils/clustalo-1.2.0"):

    # Clustal Omega (v1.2.0)
    # Multiple Sequence Alignment
    # Output : [treeid].aln alignment file and [treeid].mat distance matrix

    # output : alignment + dist matrix
    if alignment_out_path and dist_matrix_out_path:
        clustalo = ClustalOmegaCommandline(cmd=cmd_path, infile=geneSeq_file_path, outfile=alignment_out_path, distmat_full=True, distmat_out=dist_matrix_out_path,verbose=False, outfmt="clu", auto=True)
    # output : alignment
    elif alignment_out_path and not dist_matrix_out_path:
        clustalo = ClustalOmegaCommandline(cmd=cmd_path, infile=geneSeq_file_path, outfile=alignment_out_path, verbose=False, outfmt="clu", auto=True)
    # output : dist matrix
    elif not alignment_out_path and dist_matrix_out_path:
        if aligned:
            clustalo = ClustalOmegaCommandline(cmd=cmd_path, infile=geneSeq_file_path, max_hmm_iterations=-1, distmat_full=True, distmat_out=dist_matrix_out_path, verbose=False)
        else:
            clustalo = ClustalOmegaCommandline(cmd=cmd_path, infile=geneSeq_file_path, max_hmm_iterations=-1, distmat_full=True, distmat_out=dist_matrix_out_path, dealign=True, verbose=False)

    clustalo()



def phyml(geneSeq_file_path, trees_processed, treeid, cmd_path='utils/phyml'):

    input_trees = "utils/tmp/%s.newick" %treeid

    # PhyML (v20140520)
    # Calculate the log likelihood of the output tree(s) and the given gene sequence data
    with open(input_trees, "w") as newick_file:
        for tree in trees_processed:
            # Write newick for trees in a .newick file named after the treeid of the first tree
            # Convert all tree labels to gene names only
            t_tmp = tree.copy()
            leaves = t_tmp.get_leaves()
            for leaf in leaves:
                leaf.name=leaf.genes
            newick_file.write(t_tmp.write(features=[])+"\n")

    # Set everything up to run PhyML on the sequences and get the log likelihood for tree
    phyml = _Phyml.PhymlCommandline(cmd=cmd_path, input=geneSeq_file_path, input_tree=input_trees, optimize="none", bootstrap=0)
    phyml()

    # Get file extension (nexus or phylip)
    align_extension = os.path.splitext(geneSeq_file_path)[1]

    # PhyML output
    output_stats = "utils/tmp/%s%s_phyml_stats.txt" %(treeid, align_extension)
    output_tree = "utils/tmp/%s%s_phyml_tree.txt" %(treeid, align_extension)

    # Parse and fetch log-likelihood
    ll_keyword = ". Log-likelihood:"
    log_index = 0
    with open(output_stats) as search:
        for line in search:
            if ll_keyword in line:
                line = line.replace(ll_keyword, "")
                trees_processed[log_index].log_likelihood=line.strip()
                trees_processed[log_index].tree_number=log_index+1
                log_index += 1

    # Clean up tmp files
    os.remove(geneSeq_file_path)
    os.remove(input_trees)
    os.remove(output_stats)
    os.remove(output_tree)

    return trees_processed


def nexrepair(nxsfile):
    """Repair nexus file for phyML"""
    #print "REFORMATING your nexus file to phyML input ..."
    with open(nxsfile, 'r') as infile:
        with open("tmp", 'w') as outfile:
            found_matrix = False
            first_block_passed = False
            newFileContent = ""
            for line in infile:
                line = line.replace("\n", "").strip()
                lineup = line.upper()
                ignoreLine = False
                correctedLine = line
                if 'format missing=?'==line:
                    ignoreLine = True
                    #CLustal missing line found

                elif lineup == "MATRIX":
                    found_matrix = True
                    #MATRIX LINE found

                if found_matrix and not first_block_passed and line == "":
                    first_block_passed = True

                if line != "" and first_block_passed:
                    parts = line.split()
                    correctedLine = parts[-1]

                if line==";":
                    newFileContent += ";"
                elif not ignoreLine:
                    if newFileContent != "":
                        newFileContent += "\n"
                        newFileContent += correctedLine

                outfile.write(newFileContent)

    subprocess.call(['mv', "tmp", nxsfile])
    return nxsfile
