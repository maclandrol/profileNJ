import sys
from pprint import pprint

def compareFasta(ffile1,ffile2):
	fautif=''
	with open(ffile1, 'r') as FILE1, open(ffile2, 'r') as FILE2:
		fasta1 = readFasta(FILE1)
		fasta2= readFasta(FILE2)
		identiq=True
		for key in fasta1.keys():
			if fasta2[key]!=fasta1[key]:
				identiq=False
				fautif='>Diff at %s\nFile1: %s\n\nFile2: %s'%(key, fasta1[key], fasta2[key])
				break
		if(len(fasta1)== len(fasta2) and identiq):
			return True
	print fautif
	return False


#raxmlHPC -s sequenceFileName -n -d outpur.tree -m GTRCAT -b seed -p -# 100/autoMR/autoFC -f a/d D
def readFasta(infile):
	seq={}
	current=''
	for line in infile:
		if line.startswith('>'):
			current=line.strip()
			seq[current]=''
		
		elif line and current!='':
			seq[current]+=line.strip()
	return seq



if len(sys.argv)!=3:
	raise Error('Not enough input')

else:
	print compareFasta(sys.argv[1], sys.argv[2])
