#!/bin/sh 
currdir=$PWD

#for d in `ls -d $1/*/`; do

for d in $(seq 0 100); do
#	head -1 "$d$fasta"
	#seq 0 100 | parallel -j 4 -k --workdir $PWD  ./cmdexec.sh {} $d $currdir
	./cmdexec.sh $d $1 $currdir;
done
