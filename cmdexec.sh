#!/bin/sh 

fasta="$2/$1/$1.fasta"
align="$2/$1/$1.align"
clout="$2/$1/$1.fa"
outdir="$3/$2/$1"
seed=$(( $1+1 ))

rm "$2/$1"/RAxML_*
#outdir="/home/manu/Dropbox/Stage_Diro/Tree/$2$1/"
#if [ $(ls -1  $del_dir |wc -l) -ne 16 ]
#then 
#	oldraxml="$2/$1/RAxML_*"
#	rm $oldraxml
#	raxmlHPC-SSE3 -s $align -n bootstrap.align.tree -m GTRGAMMA -p $seed -x $seed -f a -# autoMR -w $outdir
#	raxmlHPC-SSE3 -s $clout -n bootstrap.fa.tree -m GTRGAMMA -p $seed -x $seed -f a -# autoMR -w $outdir

#else
#	echo "\nThis directory is already Done!!\n\n"
#fi

#clustalo -i $fasta -o $clout --outfmt=fasta --force
#RAXML best tree for 10 ML run with 
#raxmlHPC-SSE3 -s $align -n align.tree -m GTRGAMMA -p $seed -w $outdir
#raxmlHPC-SSE3 -s $clout -n fa.tree -m GTRGAMMA -p $seed -w $outdir

#Raxml rapid  Bootstrap tree with 
raxmlHPC-SSE3 -s $align -n bootstrap.align.tree -m GTRGAMMA -p $seed -x $seed -f a -# autoMR -w $outdir
#raxmlHPC-SSE3 -s $clout -n bootstrap.fa.tree -m GTRGAMMA -p $seed -x $seed -f a -# autoMR -w $outdir
