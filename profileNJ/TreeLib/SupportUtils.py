# This file is part of profileNJ
#
# Date: 02/2014

from __future__ import print_function
import subprocess
import TreeUtils
from TreeClass import TreeClass
import re
import sys
import os
import time
import glob
import numpy as np
try:
    from lxml import etree
except ImportError:
    try:
        import xml.etree.cElementTree as etree
    except ImportError:
        try:
            # Python 2.5
            import xml.etree.ElementTree as etree
        except:
            pass

def timeit(func):

    def timed(*args, **kw):
        tstart = time.time()
        result = func(*args, **kw)
        tend = time.time()
        ttime = tend - tstart
        # print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
        return ttime, result

    return timed


@timeit
def getDistMatrix(alignfile, outfile):
    cmd = "fastdist -I fasta -O phylip -e -o %s %s" % (outfile, alignfile)
    executeCMD(cmd)
    return outfile


def name_extractor(filename, ext="."):
    prefix, ext, suffix = filename.partition(ext)
    return prefix, suffix


def selectBestTree(values):
    return values.index(max(values))


@timeit
def runPolytomySolver(mltree, smap, spectree, outfile, distmat, r_option, slimit, plimit, algo):
    cmd = "python profileNJ -s %s -S %s -g %s -d %s -o %s -n -r %s --slimit=%s --plimit=%s -c %s" % (
        spectree, smap, mltree, distmat, outfile, r_option, slimit, plimit, algo)
    executeCMD(cmd)


@timeit
def runTreeFix(mltree, smap, spectree, align_ext, mltree_ext, logfile):
    cmd = "treefix -s %s -S %s -A '%s' -o %s -n %s -V 2 -l %s %s" % (
        spectree, smap, align_ext, mltree_ext, ".treefix.tree", logfile, mltree)
    executeCMD(cmd)


def getRFval(refTree_path, tree_path, unroot=False):
    refTree = TreeClass(refTree_path)
    tree = TreeClass(tree_path)
    if(unroot):
        refTree.unroot()
    rf, max_rf, c, p1, p2 = refTree.robinson_foulds(
        tree, unrooted_trees=unroot)
    return rf, max_rf


def fix_ps_out(ps_out, mltreesfile):
    nsol = 0
    dl_cost = np.inf
    with open(ps_out, 'r') as infile, open(mltreesfile, 'w+') as mlinfile:
        for line in infile:
            if(line.startswith('>')):
                dl_cost = int(line.split('=')[1].strip())
            else:
                nsol += 1
                mlinfile.write(line)
                with open("%s%s" % (ps_out, nsol), 'w') as outfile:
                    outfile.write(line)
    return nsol, dl_cost


def write_al_in_nxs_file(fastafile, outnxs="seq.nxs", al=0):

    if al:
        print()
        print('CONVERTING your fasta file to nexus format ... with "clustalw"')
        subprocess.call(['clustalw', ('-INFILE=' + fastafile),
                         ('-OUTFILE=' + outnxs), '-convert', '-OUTPUT=NEXUS'])

    else:
        print()
        print('Sequence not aligned!!')
        print('ALIGNING sequences with clustalw ...')
        clustalcmd = "clustalw -INFILE=" + fastafile + \
            " -OUTFILE=" + outnxs + " -OUTPUT=NEXUS"
        error = executeCMD(clustalcmd)
        if not error:
            print("DONE! :  sequences ALIGNED with clustalw")
        else:
            print("FAILED at clustalw execution !!!")
            return

    return nexrepair(outnxs)


def read_trees(file):
    """Read a trees file
    and yield each line
    """
    with open(file, 'r') as infile:
        for line in infile:
            if not (line.startswith('>')):
                yield line.strip()

def clustalo(geneSeq_file_path, treeid, alignment_out_path="", dist_matrix_out_path="", aligned=False, cmd_path="utils/clustalo-1.2.0"):
    from Bio.Align.Applications import ClustalOmegaCommandline
    # Clustal Omega (v1.2.0)
    # Multiple Sequence Alignment
    # Output : [treeid].aln alignment file and [treeid].mat distance matrix

    # output : alignment + dist matrix
    if alignment_out_path and dist_matrix_out_path:
        clustalo = ClustalOmegaCommandline(cmd=cmd_path, infile=geneSeq_file_path, outfile=alignment_out_path,
                                           distmat_full=True, distmat_out=dist_matrix_out_path, verbose=False, outfmt="clu", auto=True)
    # output : alignment
    elif alignment_out_path and not dist_matrix_out_path:
        clustalo = ClustalOmegaCommandline(
            cmd=cmd_path, infile=geneSeq_file_path, outfile=alignment_out_path, verbose=False, outfmt="clu", auto=True)
    # output : dist matrix
    elif not alignment_out_path and dist_matrix_out_path:
        if aligned:
            clustalo = ClustalOmegaCommandline(
                cmd=cmd_path, infile=geneSeq_file_path, max_hmm_iterations=-1, distmat_full=True, distmat_out=dist_matrix_out_path, verbose=False)
        else:
            clustalo = ClustalOmegaCommandline(cmd=cmd_path, infile=geneSeq_file_path, max_hmm_iterations=-1,
                                               distmat_full=True, distmat_out=dist_matrix_out_path, dealign=True, verbose=False)

    clustalo()


def nexrepair(nxsfile):
    """Repair nexus file for phyML"""
    # print "REFORMATING your nexus file to phyML input ..."
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
                if 'format missing=?' == line:
                    ignoreLine = True
                    # CLustal missing line found

                elif lineup == "MATRIX":
                    found_matrix = True
                    # MATRIX LINE found

                if found_matrix and not first_block_passed and line == "":
                    first_block_passed = True

                if line != "" and first_block_passed:
                    parts = line.split()
                    correctedLine = parts[-1]

                if line == ";":
                    newFileContent += ";"
                elif not ignoreLine:
                    if newFileContent != "":
                        newFileContent += "\n"
                        newFileContent += correctedLine

                outfile.write(newFileContent)

    subprocess.call(['mv', "tmp", nxsfile])
    return nxsfile


def executeCMD(cmd, dispout=False):
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    # print "STDERR\n---------\n", err
    if dispout:
        print("\nSTDOUT\n---------\n", out)
    return err

def retrieveDupAndLostCost(treefile, streefile, smap, sep=None, pos='prefix'):

    genetree = TreeClass(treefile)
    specietree = TreeClass(streefile)
    regexmap = {}
    speciemap = {}
    with open(smap, 'rU') if isinstance(smap, basestring) else smap as INPUT:
        for line in INPUT:
            g, s = line.strip().split()
            g_regex = re.compile(g.replace('*', '.*'))
            regexmap[g_regex] = s

    for leaf in genetree:
        for key, value in regexmap.iteritems():
            if key.match(leaf.name):
                speciemap[leaf.name] = value

    genetree.set_species(speciesMap=speciemap, sep=sep, pos=pos)
    lcamap = TreeUtils.lcaMapping(genetree, specietree)
    TreeUtils.reconcile(genetree, lcaMap=lcamap, lost="yes")
    # print genetree.get_ascii(show_internal=True, attributes=['name', 'type'])
    return TreeUtils.computeDLScore(genetree)


def consel(inputfile, type, in_ext='.txt', sort=9):
    makermtcmd = "makermt --%s %s" % (type, inputfile)
    basename, out = name_extractor(inputfile, ext=in_ext)
    conselcmd = "consel %s" % basename
    catpvcmd = "catpv %s > %s-pv.txt" % (basename, basename)
    if sort:
        catpvcmd += " -s %s" % sort
    executeCMD(makermtcmd)
    executeCMD(conselcmd)
    executeCMD(catpvcmd)
    conselOut = basename + "-pv.txt"
    return parseConselOutfile(conselOut, sort)


def parseConselOutfile(file, sort):
    title = []
    content = []
    first = True
    with open(file, 'r') as CONSELOUT:
        for line in CONSELOUT:
            if(line.startswith('#')):
                if first and 'reading' not in line and ("%s+" % sort not in line):
                    title.extend(
                        line.replace('#', '').replace('|', '').strip().split())
                    first = False
                elif not first:
                    values = line.replace('#', '').replace(
                        '|', '').strip().split()
                    content.append([float(x) for x in values])
    return dict(zip(title, zip(*content)))


def runRaxmlPval(basedir, alignfile, narbres, out=None, listfile=[], sort=None):
    # raxmlHPC-SSE3 -f g -z polytomysolver.tree -s 341.phy -m GTRGAMMA -n TEST
    if not out:
        out = "%s/all.tree" % basedir
        with open(out, 'w') as fout:
            for f in listfile:
                print(open(f, 'r').read().strip().replace(
                    '\s', '').replace('\n', ''), file=fout)

    # cleaning for raxml
    for f in glob.glob("%s/*trees" % basedir):
        os.remove(f)

    cmd = "raxmlHPC-SSE3 -f g -z %s -s %s -m GTRGAMMA -n trees -w %s" % (
        out, alignfile, os.path.abspath(basedir))
    rx_time, rst = runRAxML_likelihood(cmd)
    os.remove(out)
    infofile = glob.glob("%s/*info.trees" % basedir)[0]
    consel_input = glob.glob("%s/*perSiteLLs.trees" % basedir)[0]
    likelihoods = extractRAXMLikelihood(infofile, narbres)
    consel_output = {}

    title = 'rank item    obs     au     np      bp     kh     sh    wkh    wsh'
    if(narbres) > 1:
        consel_output = consel(consel_input, 'puzzle', '.trees', sort)
    else:
        consel_output = dict(zip(title.split(), 10 * ['N/A']))

    consel_output['likelihood'] = likelihoods
    return rx_time, consel_output


def calculate_likelihood(profileNJ_file, alignfile, basedir=os.getcwd()):
    cmd = "raxmlHPC-SSE3 -f G -z %s -s %s -m GTRGAMMA -n trees -w %s" % (
        profileNJ_file, alignfile, os.path.abspath(basedir))
    rx_time, rst = runRAxML_likelihood(cmd)
    infofile = glob.glob("%s/*info.trees" % basedir)[0]
    consel_input = glob.glob("%s/*perSiteLLs.trees" % basedir)[0]
    likelihoods = extractRAXMLikelihood(infofile, narbres)
    consel_output = {}

    title = 'rank item    obs     au     np      bp     kh     sh    wkh    wsh'
    if(narbres) > 1:
        consel_output = consel(consel_input, 'puzzle', '.trees', sort)
    else:
        consel_output = dict(zip(title.split(), 10 * ['N/A']))

    consel_output['likelihood'] = likelihoods
    return rx_time, consel_output


@timeit
def runRAxML_likelihood(cmd):
    executeCMD(cmd)


def extractRAXMLikelihood(filename, n):
    likelihoods = []
    with open(filename) as IN:
        patern = re.compile('Tree [0-9]+: ')
        for line in reversed(IN.readlines()):
            if (patern.match(line)):
                likelihoods.append(float(line.split(':')[1].strip()))
                n -= 1
            if(n <= 0):
                break

    return list(reversed(likelihoods))
