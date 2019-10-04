#!/bin/bash

#example: bash loop_aino_merge-bedgraphs.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT/

#save inputdirectories into variables                                                                                                                                       
alignment_dir=$1
mem=2gb    

mkdir -p $alignment_dir/log/merge-bedgraphs 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $alignment_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_ ]     
    then
        echo $dir
        qsub -N merge-bedgraphs_$name -l mem=$mem,walltime=4:00:00 -o $alignment_dir/log/merge-bedgraphs/merge-bedgraphs_${name}.log -e $alignment_dir/log/merge-bedgraphs/merge-bedgraphs_${name}.err -v xdir=${dir},filePlus=$dir/${name}_5p.dz-collapsed.plus.bedgraph,fileMinus=$dir/${name}_5p.dz-collapsed.minus.bedgraph,outfile=$dir/${name}_5p.dz-collapsed.bedgraph aino_merge-bedgraphs.sh
    fi
done
