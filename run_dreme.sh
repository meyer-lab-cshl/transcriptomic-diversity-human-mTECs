#!/bin/bash

#example: bash run_dreme.sh /home/stroemic/hiwi_16/analysis/gene_lists/sequences/aire_dep_tras_mTECs_sequences.fa /home/stroemic/hiwi_16/analysis/gene_lists/sequences/other_tras_mTECs_sequences.fa /home/stroemic/hiwi_16/analysis/meme/dreme/aire_vs_others/ 3

#save inputdirectories into variables                                                                                                                             
file=$1                                                                                                                                                
background=$2
outdir=$3
no_motifs=$4
mem=2gb
                                                                                                                                                          
mkdir -p $outdir                                                                                                                                      
mkdir -p $outdir/log                                                                                                                                      

logdir=$outdir/log

name=`echo $file | perl -pe 's/.*\/(.*)_sequences.fa/$1/'`  
                                  
if [ ! -d $outdir/dreme.html ]
then                                                                                                                                                
   echo $name
   qsub -N dreme_$name -l mem=$mem,walltime=4:00:00 -o $logdir/dreme_${name}.log -e $logdir/dreme_${name}.err -v file=$file,background=$background,outdir=$outdir,no_motifs=$no_motifs dreme.sh
fi                                                                                                                                                 

