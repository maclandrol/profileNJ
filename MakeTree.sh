#!/bin/zsh 
currdir=$PWD

#for d in `ls -d $1/*/`; do
biggest=( $(for dossier in $(seq 0 999); do echo "`grep '>' "$1/$dossier/$dossier.fasta"|wc -l` $dossier"; done | sort -g -s -r -k 1 | awk '{print $2}') )
echo $biggest > "$1/order.txt"
for d in $(seq 40 50); do
#	head -1 "$d$fasta"
	#seq 0 100 | parallel -j 4 -k --workdir $PWD  ./cmdexec.sh {} $d $currdir
	./cmdexec.sh ${biggest[$d]} $1 $currdir;
done