#!/bin/bash

#example: bash loop_leonie_properBedgraph.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT

#save inputdirectories into variables                                                                                                                                         
bedgraph_dir=$1
mem=5gb    

mkdir -p $bedgraph_dir/log/properBedgraph 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $bedgraph_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        
    
    if [ ! -f $dir/${name}_5p.dz-collapsed.upload.bedgraph ]     
    then

        qsub -N properBedgraph_$name -l mem=$mem,walltime=4:00:00 -o $bedgraph_dir/log/properBedgraph/properBedgraph_${name}.log -e $bedgraph_dir/log/properBedgraph/properBedgraph_${name}.err -v xdir=$dir,name=$name,bedgraph=$dir/${name}_5p.dz-collapsed.bedgraph,outfile=$dir/${name}_5p.dz-collapsed.converted.bedgraph leonie_properBedgraph.sh
    fi
done
