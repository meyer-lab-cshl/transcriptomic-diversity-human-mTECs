#!/bin/bash

#example: bash loop_meme.sh /home/stroemic/hiwi_16/analysis/gene_lists/sequences/ /home/stroemic/hiwi_16/analysis/meme/ oops(or anr/zoops) 2

#save inputdirectories into variables                                                                                                                             
seqdir=$1                                                                                                                                                
outdir=$2
mode=$3
no_motifs=$4
mem=2gb
                                                                                                                                                          
mkdir -p $outdir                                                                                                                                      
mkdir -p $outdir/log                                                                                                                                      

logdir=$outdir/log

#iterate over all files in fastqdir, extract file identifier, run fatsqc on file and write output into outputdir, directly extract files from zip archiven
for file in $seqdir/*.fa; do                                                                                                                    

    if [[ $file =~ hk.*.fa ]] || [[ $file =~ other.*.fa ]] #make sure to exclude hk and other_tras classes, they are to big 
    then
        continue
    fi

    name=`echo $file | perl -pe 's/.*\/(.*)_sequences.fa/$1/'`  
                                  
    if [ ! -d $outdir/${name}/bla ]                                                                                                                   
    then                                                                                                                                                
        echo $name
        qsub -N meme_$name -l mem=$mem,walltime=6:00:00 -o $logdir/meme_${name}.log -e $logdir/meme_${name}.err -v outdir=$outdir/${name},file=$file,mode=$mode,no_motifs=$no_motifs meme.sh
    fi                                                                                                                                                 
done 
