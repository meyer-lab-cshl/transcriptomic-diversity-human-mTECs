#!/bin/bash

#example: bash loop_unique_reads_sam.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
mem=5gb    

mkdir -p $sam_dir/log/unique_reads 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $sam_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_Aligned.out.unique.sorted.bam ]     
    then
        qsub -N unique_bam_$name -l mem=$mem,walltime=2:00:00 -o $sam_dir/log/unique_reads/unique_reads_${name}.log -e $sam_dir/log/unique_reads/unique_reads_${name}.err -v filedir=$dir,name=$name unique_reads_sam.sh
    fi
done
