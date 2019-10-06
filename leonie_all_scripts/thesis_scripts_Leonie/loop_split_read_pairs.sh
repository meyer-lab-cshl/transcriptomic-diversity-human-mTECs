#!/bin/bash

#example: bash loop_split_read_pairs.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
mem=2gb    

mkdir -p $sam_dir/log/split_read_pairs 

#iterate over all subdirs in inputdir, extract identifier, call script which uses samtools with flags to extract single reads from pairs
for dir in $sam_dir/C*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_Aligned.out.unique.sorted.mbc._1.sam ]     
    then
        qsub -N split_read_pairs_$name -l mem=$mem,walltime=2:00:00 -o $sam_dir/log/split_read_pairs/split_read_pairs_${name}.log -e $sam_dir/log/split_read_pairs/split_read_pairs_${name}.err -v filedir=$dir,name=$name split_read_pairs.sh
    fi
done
