from PolytomySolver.Multipolysolver import *
from TreeLib import *
from TreeLib.SupportUtils import *
import os, shutil, glob

#insert notung in our comparision
#what to compare :
#For the same alpha parameters in 
#number of solution found (same)
# la resolution de polytomie de notung essaye de maintenir le DL le plus bas pour tous l'Arbre alors que pour nous c'est pour la polytomie uniquement, en se basant sur le fait
#que la resolution independante de polytomie donne l'Arbre avec e score mini

@timeit
def notungResolve(basedir, genetree, specietree, seuil, task="rearrange"):
	"""Execute Notung"""
	notung="~/Programmes/Notung-2.6/Notung-2.6.jar"
	output=[]
	rootcmd="java -jar %s -g %s -s %s --root --costdup 1  --outputdir %s"%(notung, genetree, specietree, basedir)
	executeCMD(rootcmd)
	outputfile="%s.rooting.0"%genetree
	
	if(task=="rearrange"):
		rearrangecmd= "java -jar %s -g %s -s %s --rearrange --threshold %s --treeoutput newick  --costdup 1 --nolosses --maxtrees 100 --speciestag prefix  --outputdir %s"%(notung, outputfile, specietree, seuil, basedir)
		executeCMD(rearrangecmd)
		output= glob.glob("%s/*rearrange*"%basedir)

	elif(task=="resolve"):
		resolvecmd= "java -jar %s -g %s -s %s --resolve --treeoutput newick  --costdup 1 --nolosses --maxtrees 100 --speciestag prefix  --outputdir %s"%(notung, outputfile, specietree, basedir)
		output= glob.glob("%s/*resolve*"%basedir)

	return output

		
def headerWrite(stream, header):
	line ="\t".join(header)+"\n"
	if line not in stream.readline():
		stream.write(line)


def selectBestTree(values):
	return values.index(max(values))

def runTest(outfile, basedir, seq, smap, specietree, seuil):
	deletethis = []
	print "Cleaning process ..."
	deletethis.extend(glob.glob("%s/*rooting*"%basedir))
	deletethis.extend(glob.glob("%s/*perSiteLLs.trees"%basedir))
	for file in set(deletethis):
		os.remove(file)
	
	inputtree= os.path.join(basedir, "%s.align.bootstrap.tree"%(seq))
	truetree=os.path.join(basedir, "%s.tree"%(seq))
	alignfile=os.path.join(basedir, "%s.align"%seq)
	header=['tree', 'Notung_dlc','Notung_nad','Notung_dup','Notung_loss', 'Notung_rf','Notung_maxrf', 'time', 'ML_time', 'notung_nsol', 'bestlkl']
	for s in seuil:
		with open("%s%s_notung.csv"%(outfile, s), 'a+') as OUT:
			headerWrite(OUT, header)
			nt_time, treesfile= notungResolve(basedir, inputtree, specietree, s, task="rearrange")
			mltime, consel_out= runRaxmlPval(basedir, alignfile, len(treesfile), out=None, listfile=treesfile, sort=9)
			bestposition= int(consel_out['item'][0])-1 if len(treesfile)>1 else selectBestTree(consel_out['likelihood'])
			besttree=treesfile[bestposition]
			notung_rf, notung_maxrf=getRFval(truetree, besttree)
			notung_nad, notung_dup, notung_loss= retrieveDupAndLostCost(besttree, specietree, smap, sep='_', pos='prefix')
			dlc=notung_nad+notung_dup+notung_loss
			line=[seq, dlc, notung_nad, notung_dup, notung_loss,notung_rf,notung_maxrf, nt_time, mltime, len(treesfile), consel_out['likelihood'][bestposition]]
			OUT.write("\t".join([str(val) for val in line])+"\n")

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Simulation script, run polytomySolver,TreeFix and compare the output. The input file name should respect TreeFix standard')
	parser.add_argument('-s', '--species', dest='specietree', help="Name of the file containing the species newick tree.",required=True)
	parser.add_argument('-S', '--smap', dest='smap', help="Gene to species map. Use the standard format.",required=True)
	parser.add_argument('-g', '--mltree', dest='mltree', help="Name of the file containing the gene newick tree.")
	parser.add_argument('-w', '--workdir', dest='workdir', default="./", help="Working directory. The directory which contains each simulation")
	parser.add_argument('--seqlist', dest='seqlist', nargs=2, help="Choose the file listing your input and the number of input to process.",required=True)
	parser.add_argument('--seuil', dest='seuil', type=int, nargs='+', help="Support contraction threshold for polytomySolver")
	parser.add_argument('--out', dest='out', default="output.csv", help="Output file of the analysis.")
	args = parser.parse_args()

	outfile=os.path.join(args.workdir, args.out)
	instream=open(args.seqlist[0], 'r').read().strip().split()
	seq_list=[int(x) for x in instream][:int(args.seqlist[1])]
	def_tree="RAxML_bipartitions.bootstrap.align.tree"

	for seq in seq_list:
		basedir=os.path.join(args.workdir, str(seq))
		if(os.path.exists(os.path.join(basedir, def_tree))):
			runTest(outfile, basedir, seq, args.smap, args.specietree, args.seuil)
	
