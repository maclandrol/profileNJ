""" 
This script is a pipeline to use with the TreeClass tools.
........................... 
Missing doc!!!!!!!!!!
"""
import TreeUtils
import subprocess
from TreeClass import TreeClass
import re
import numpy as np

phymllk='phyml_lk.txt'
phymltree='_phyml_tree.txt'
phymlmodel='_phyml_stats.txt'
phymlbootstrap='_phyml_boot_trees.txt'
phymlbootstrapmodel='_phyml_boot_stats.txt'
phymltrees='_phyml_trees.txt'

def executePipe(tree, nxsfile=None, fasta=None, al=0, type=None, treefile=None):
	
	n=[]
	for leaf in tree:
		if(len(n)<7):
			n.append(leaf.name)
	
	tree.prune(n)

	if(treefile is not None):
		tree=TreeClass(treefile)
	else:
		try:
			tree.write(format=0, outfile="tree.nw");
			treefile="tree.nw"
		except Exception as e:
			print e
			print "Can't write tree to 'tree.nw'"

	if not isinstance(tree, TreeClass):
		raise ValueError ("You sould use a TreeNode instance")

	if(nxsfile is None):
		if fasta is None:
			print 
			print "WRITING your sequence into a fasta file"
			tree.writeSeqToFasta(comment=0)
			fasta="seq.fasta"
	
		nxsfile=write_al_in_nxs_file(fasta, al=al)

	executePhyML(nxsfile, treefile)


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
			print "SUCCESSFULLY DONE! :  sequences ALIGNED with clustalw"
		else:
			print "FAILED at clustalw execution !!!"
			return


	return nexrepair(outnxs)


def executePhyML(seqfile, treefile, quiet=0):
	print 
	print "EXECUTING phyML ..."
	phymlcmd= "phyml -i " + seqfile +" -u " + treefile + " -n 1 -o l --print_site_lnl --no_memory_check"
	if(quiet):
		phymlcmd=phymlcmd+ " --quiet"
	print phymlcmd
	error=executeCMD(phymlcmd)
	if not error:
		print "SUCCESSFULLY DONE! : phyML EXECUTED"
	else:
		print "FAILED at phyML execution !!!"
		return


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


"""I'm ashamed of this, truly ashamed"""
def nexrepair(nxsfile):
	print "REFORMATING your nexus file to phyML input ..."
	with open(nxsfile, 'r') as infile, open("tmp", 'w') as outfile:
		inMatrix = False
		firstBlockPassed = False
		newFileContent = ""
		for line in infile:
			line = line.replace("\n", "")
		  	lineup = line.upper()
			ignoreLine = False
			correctedLine = line 
			if line == 'format missing=?':
				ignoreLine = True
				print "CLustal missing line found"
			
			elif lineup == "MATRIX":
				inMatrix = True
				print "MATRIX LINE found"
		    
			if inMatrix and not firstBlockPassed and line == "":
				firstBlockPassed = True
		    
			if line != "" and firstBlockPassed:
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


def capitalize(line):
	return "".join([line[0].upper(), line[1:]])
	

def customTreeCompare(original_t, corrected_t, t):

	#Leaves remaining test and original binary node test
	ct_leaves=[]
	t_leaves=[]
	t_binary=[]
	success=[]
	ct_binary=[]
	for node in original_t.traverse("levelorder"):
		desc_name=set(node.get_leaf_names())
		ct_parent= corrected_t.get_common_ancestor(desc_name)
		t_parent= t.get_common_ancestor(desc_name)
		ctl=set(ct_parent.get_leaf_names())
		tl=set(t_parent.get_leaf_names())
		ct_leaves.append(ctl.difference(desc_name))
		t_leaves.append(tl.difference(desc_name))
		if(node.is_binary() and not node.has_polytomies()):
			ct_binary.append(ct_parent.robinson_foulds(node)[0:3])
			t_binary.append(t_parent.robinson_foulds(node)[0:3])
		success.append(len(tl.difference(ctl))<1)



	ct_success=filter(None, map(lambda x:len(x)<1,ct_leaves))
	t_success=filter(None, map(lambda x:len(x)<1,t_leaves))

	print "\nCorrected Tree binary list rf_fould\n"
	print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_binary))
	print "\nTree binary list rf_fould\n"
	print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), t_binary))

	if(len(ct_success)==len(ct_leaves)):
		print "**Leave remaining success for corrected tree"
		#print "\n".join([str(h) for h in t_success])
	else:
		print "**Corrected tree doesn't follow patern"
		#print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_leaves))


	if(len(t_success)==len(t_leaves)):
		print "**Leave remaining success for tree"
		#print "\n".join([str(h) for h in t_success])
	else:
		print "**Tree doesn't follow patern"
		#print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), t_leaves))

	print "**Compatibility test between tree: ", all(success)
	
if __name__ == '__main__':
	ensemblTree=TreeUtils.fetch_ensembl_genetree_by_id(treeID="ENSGT00390000003602", output="phyloxml", aligned=0, sequence="cdna")
	executePipe(ensemblTree, al=0, type="cdna")
