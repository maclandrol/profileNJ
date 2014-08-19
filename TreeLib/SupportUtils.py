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
        return result, ttime

    return timed


def methodCompare(mltree, smap, specietree, alignfile, gtree, seuil, mltree_ext, r_option, slimit, plimit, correctPhylo):
    w_dir=os.path.dirname(mltree)
    basename, align_ext=name_extractor(os.path.basename(alignfile), ext=".")
    distmat= getDistMatrix(alignfile,os.path.join(w_dir, basename+".dist"))
    logfile= os.path.join(w_dir, basename+".treefix.log")
    ps_out=os.path.join(w_dir, basename+".polytomysolver.tree")
    tf_out=os.path.join(w_dir, basename+".treefix.tree")

    runPolytomySolver(gtree, smap, specietree, ps_out, distmat, r_option, slimit, plimit)
    runTreeFix(mltree, smap, specietree, "."+align_ext, mltree_ext, logfile)
    fix_ps_out(ps_out)
    # compute pval and Dlnl for true tree using RAxML tree to optimize
    tf_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s %s" %("."+align_ext, ".treefix.tree", ".tf.ml", correctPhylo)
    executeCMD(tf_cmp_lk)
    ps_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s %s" %("."+align_ext, ".polytomysolver.tree",".ps.ml", correctPhylo)
    executeCMD(ps_cmp_lk)

    tt_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s %s" %("."+align_ext, ".tree",".tt.ml", correctPhylo)
    executeCMD(tt_cmp_lk)

    rx_cmp_lk= "treefix_compute --type likelihood -m treefix.models.raxmlmodel.RAxMLModel -A %s -U %s -n %s %s" %("."+align_ext,mltree_ext ,".rx.ml", correctPhylo)
    executeCMD(rx_cmp_lk)

    # compute dl cost
    raxml_cmp_dl="treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, mltree_ext, ".raxml.output", mltree)
    executeCMD(raxml_cmp_dl)

    tf_cmp_dl="treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".treefix.tree", ".treefix.output", tf_out)
    executeCMD(tf_cmp_dl)

    ps_cmp_dl="treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".polytomysolver.tree", ".polytomysolver.output", ps_out)
    executeCMD(ps_cmp_dl)

    tt_cmp_dl="treefix_compute --type cost -m treefix.models.duplossmodel.DupLossModel  -s %s -S %s -o %s -n %s %s"%(specietree, smap, ".tree", ".true.output", correctPhylo)
    executeCMD(tt_cmp_dl)

def getDistMatrix(alignfile, outfile):
    cmd="fastdist -I fasta -O phylip -e -o %s %s"%(outfile, alignfile)
    executeCMD(cmd)
    return outfile


def name_extractor(filename, ext="."):
    prefix, ext, suffix = filename.partition(ext)
    return prefix, suffix


#raxmlHPC-SSE3 -f I -m GTRGAMMA -t 0.bootstrap.align.tree -n broot
#
@timeit
def runPolytomySolver(mltree, smap, spectree, outfile, distmat, r_option, slimit, plimit):
    cmd="python PolytomySolver.py -s %s -S %s -g %s -d %s -o %s -n -r %s --slimit=%s --plimit=%s"%(spectree, smap, mltree, distmat, outfile, r_option, slimit, plimit)
    executeCMD(cmd)


@timeit
def runTreeFix(mltree, smap, spectree, align_ext, mltree_ext, logfile):
    cmd="treefix -s %s -S %s -A '%s' -o %s -n %s -V 2 -l %s %s"%(spectree, smap, align_ext, mltree_ext, ".treefix.tree", logfile, mltree)
    executeCMD(cmd)


def fix_ps_out(ps_out):
    with open(ps_out, 'r') as infile, open("tmp", 'w') as outfile:
        keepline=''
        for line in infile:
            if(line.startswith('>')):
                print line
            else:
                keepline=line
        outfile.write(keepline)
    try:
        os.rename("tmp", ps_out)
    except Exception, e:
        pass


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
