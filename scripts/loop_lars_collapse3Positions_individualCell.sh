#!/bin/bash

#example: bash loop_lars_collapse3Positions_individualCell.sh /home/stroemic/hiwi_16/data/Summary_counts/ /home/stroemic/hiwi_16/analysis/lars_perl/ /home/stroemic/hiwi_16/analysis/lars_perl/all_consensus_sites.csv

#save inputdirectories into variables                                                                                                                                         
sc_dir=$1
out_dir=$2
isoform_deffile=$3
mem=10gb    

mkdir -p $out_dir/log/individual_samples 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for file in $sc_dir/*.summary.counts ; do                          

    name=`echo $file | perl -pe 's/.*\/(.*)\.summary.counts/$1/'`   
    echo $name        

    if [ ! -f $out_dir/${name}.genes.csv ]     
    then

        qsub -N collapse_individual_$name -l mem=$mem,walltime=2:00:00 -o $out_dir/log/individual_samples/individual_samples_${name}.log -e $out_dir/log/individual_samples/individual_samples_${name}.err -v sc_file=$file,isoform_deffile=$isoform_deffile,out_file_isoforms=$out_dir/${name}.isoforms.csv,out_file_genes=$out_dir/${name}.genes.csv lars_collapse3Positions_individualCell.sh
    fi
done
