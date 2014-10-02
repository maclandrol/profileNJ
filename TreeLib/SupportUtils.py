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


def timeit(func):

    def timed(*args, **kw):
        tstart = time.time()
        result = func(*args, **kw)
        tend = time.time()
        ttime= tend-tstart

        print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
        return ttime, result

    return timed


def getDistMatrix(alignfile, outfile):
    cmd="fastdist -I fasta -O phylip -e -o %s %s"%(outfile, alignfile)
    executeCMD(cmd)
    return outfile


def name_extractor(filename, ext="."):
    prefix, ext, suffix = filename.partition(ext)
    return prefix, suffix

@timeit
def phymlikelihood(sequence, align_ext, treepath, n_sol, s_model='HKY85'):
    
    sequence=convert_to_phylip(sequence, 'fasta', 'align')
    likelihoods=[]
    output_stats = "%s_phyml_stats.txt" %sequence
    output_tree = "%s_phyml_tree.txt" %sequence
    ll_keyword = ". Log-likelihood:"
 
    for n in xrange(n_sol):
        phyml = _Phyml.PhymlCommandline(cmd="phyml", input=sequence, input_tree="%s%s"%(treepath,(n+1)) , optimize="none", bootstrap=0, model=s_model)
        phyml()

        with open(output_stats) as search:
            for line in search:
                if ll_keyword in line:
                    line = line.replace(ll_keyword, "")
                    likelihoods.append(float(line))

    return likelihoods

def selectBestTree(likelihoods):
    return likelihoods.index(max(likelihoods))


@timeit
def runPolytomySolver(mltree, smap, spectree, outfile, distmat, r_option, slimit, plimit, algo):
    cmd="python PolytomySolver.py -s %s -S %s -g %s -d %s -o %s -n -r %s --slimit=%s --plimit=%s -c %s"%(spectree, smap, mltree, distmat, outfile, r_option, slimit, plimit, algo)
    executeCMD(cmd)


@timeit
def runTreeFix(mltree, smap, spectree, align_ext, mltree_ext, logfile):
    cmd="treefix -s %s -S %s -A '%s' -o %s -n %s -V 2 -l %s %s"%(spectree, smap, align_ext, mltree_ext, ".treefix.tree", logfile, mltree)
    executeCMD(cmd)


def getRFval(refTree_path, tree_path, unroot=False):
    refTree=TreeClass(refTree_path)
    tree=TreeClass(tree_path)
    rf, max_rf, c, p1, p2= refTree.robinson_foulds(tree, unrooted_trees=unroot)
    return rf, max_rf


def fix_ps_out(ps_out):
    nsol=0
    dl_cost=np.inf
    with open(ps_out, 'r') as infile:
        for line in infile:
            if(line.startswith('>')):
                dl_cost=int(line.split('=')[1].strip())
            else:
                nsol+=1
                with open("%s%s"%(ps_out, nsol), 'w') as outfile:
                    outfile.write(line)
    return nsol, dl_cost


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
    phyml = _Phyml.PhymlCommandline(cmd=cmd_path, input=geneSeq_file_path, input_tree=input_trees, optimize="none", bootstrap=-1)
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


def executeCMD(cmd):
    print "\n", cmd
    p=subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    print "STDERR\n---------\n" , err
    print "\nSTDOUT\n---------\n", out
    return err
