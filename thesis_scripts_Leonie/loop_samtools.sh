#!/bin/bash

#example: bash loop_samtools.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT plus

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
strand=$2
mem=5gb    

mkdir -p $sam_dir/log/samtools 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $sam_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_Aligned.out.unique.sorted.mbc.collapsed_1.${strand}.sorted.bam ]     
    then
        qsub -N samtools_$name -l mem=$mem,walltime=18:00:00 -o $sam_dir/log/samtools/samtools_${name}.log -e $sam_dir/log/samtools/samtools_${name}.err -v filedir=$dir,name=$name,strand=$strand samtools.sh
    fi
done
