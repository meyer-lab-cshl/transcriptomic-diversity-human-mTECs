#!/bin/bash

#example: bash loop_fastqc_on_fastq.sh /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/ /home/stroemic/hiwi_16/analysis/FastQC/FastQC_5Seq_processed_rerun

#save inputdirectories into variables                                                                                                                             
fastqdir=$1                                                                                                                                                
fastqcresultsdir=$2                                                                                                                                              
mem=1gb
                                                                                                                                                          
mkdir -p $fastqcresultsdir                                                                                                                                      
mkdir -p $fastqcresultsdir/log                                                                                                                                      

logdir=$fastqcresultsdir/log

#iterate over all files in fastqdir, extract file identifier, run fatsqc on file and write output into outputdir, directly extract files from zip archiven
for file in $fastqdir/*.fastq; do                                                                                                                    

    #if [[ $file =~ .*_2.fastq ]]  #make sure to exclude one  pair end read
    #then
    #    continue
    #fi

    name=`echo $file | perl -pe 's/.*\/(.*)\.fastq/$1/'`  
                                  
    if [ ! -d $fastqcresultsdir/${name}_fastqc ]                                                                                                                   
    then                                                                                                                                                
        echo $name
    qsub -N fastqc_processed_$name -l mem=$mem,walltime=3:00:00 -o $logdir/fastqc_processed_rerun_${name}.log -e $logdir/fastqc_processed_rerun_${name}.err -v outdir=$fastqcresultsdir,file=$file fastqc_on_fastq.sh
    fi                                                                                                                                                 
done 
