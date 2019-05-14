#!/bin/bash

#example: bash loop_collapse_5pos.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT/ plus

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
strand=$2
mem=5gb    

mkdir -p $sam_dir/log/collapse_5pos 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $sam_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_5p.dz-collapsed.${strand}.bedgraph ]     
    then

        qsub -N collapse_5pos_$name -l mem=$mem,walltime=2:00:00 -o $sam_dir/log/collapse_5pos/collapse_5pos_${name}.log -e $sam_dir/log/collapse_5pos/collapse_5pos_${name}.err -v bamfile=$dir/${name}_Aligned.out.unique.${strand}.sorted.bam,outfile=$dir/${name}_5p.dz-collapsed.${strand}.bedgraph collapse_5pos.sh
    fi
done
