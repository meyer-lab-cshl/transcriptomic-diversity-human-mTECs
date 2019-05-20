#!/bin/bash

#example: bash loop_aino_divideByStrand.sh /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT/

#save inputdirectories into variables                                                                                                                                         
sam_dir=$1
mem=2gb    

mkdir -p $sam_dir/log/divideByStrand 

#iterate over all subdirs in inputdir, extract identifier, use samtools view with flag to extract unique identifiers
for dir in $sam_dir/*; do                          
    if [[ ! $dir =~ C.* ]] 
    then
        continue
    fi

    name=`echo $dir | perl -pe 's/.*\/(.*)/$1/'`   
    echo $name        

    if [ ! -f $dir/${name}_Aligned.out.unique.plus.sam ]     
    then

        qsub -N divideByStrand_$name -l mem=$mem,walltime=2:00:00 -o $sam_dir/log/divideByStrand/divideByStrand_${name}.log -e $sam_dir/log/divideByStrand/divideByStrand_${name}.err -v samfile=$dir/${name}_Aligned.out.unique.sam,ofilePlus=$dir/${name}_Aligned.out.unique.plus.sam,ofileMinus=$dir/${name}_Aligned.out.unique.minus.sam aino_divideByStrand.sh
    fi
done
