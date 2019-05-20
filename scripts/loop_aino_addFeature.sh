#!/bin/bash

#example: bash loop_aino_addFeature.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT/ /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/mol_bc

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
molbc_dir=$2
mem=45gb    

mkdir -p $sam_dir/log/addFeature 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $sam_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_Aligned.out.unique.sorted.mbc.sam ]     
    then

        qsub -N addFeature_$name -l mem=$mem,walltime=4:00:00 -o $sam_dir/log/addFeature/addFeature_${name}.log -e $sam_dir/log/addFeature/addFeature_${name}.err -v samfile=$dir/${name}_Aligned.out.unique.sorted.sam,molbctsv=$molbc_dir/${name}.molbc.tsv,outfile=$dir/${name}_Aligned.out.unique.sorted.mbc.sam aino_addFeature.sh
    fi
done
