#!/bin/bash

#example: bash loop_aino_5Seq_rename_collapse.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT/ 1 

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
no=$2
mem=5gb    

mkdir -p $sam_dir/log/collapsing_$no 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $sam_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name 
    echo $no

    if [ ! -f $dir/${name}_Aligned.out.unique.sorted.mbc.collapsed_${no}.sam ]     
    then
        mkdir -p $dir/collapsing_${no} 
        qsub -N collapsing_$name -l mem=$mem,walltime=4:00:00 -o $sam_dir/log/collapsing_${no}/collapsing_${name}.log -e $sam_dir/log/collapsing_${no}/collapsing_${name}.err -v samfile=$dir/${name}_Aligned.out.unique.sorted.mbc._${no}.sam,collapsed_file=$dir/${name}_Aligned.out.unique.sorted.mbc.collapsed_${no}.sam,mc_count_file=$dir/collapsing_${no}/${name}_collapsed_countsPerMol.tsv,pos_count_file=$dir/collapsing_${no}/${name}_collapsed_countsPerPos.tsv,collapse_info=$dir/collapsing_${no}/${name}_collapsed.info aino_5Seq_rename_collapse.sh
    fi
done
