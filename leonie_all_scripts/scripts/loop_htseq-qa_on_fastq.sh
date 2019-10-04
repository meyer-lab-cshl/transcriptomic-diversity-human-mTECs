#!/bin/bash

#example: bash loop_htseq-qa_on_fastq.sh /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed  /home/stroemic/hiwi_16/analysis/htseq-qa/5Seq_processed
#save inputdirectories into variables                                                                                                                             
fastqdir=$1                                                                                                                                                
resultsdir=$2                                                                                                                                              
mem=1gb
                                                                                                                                                          
mkdir -p $resultsdir                                                                                                                                      
mkdir -p $resultsdir/log                                                                                                                                      

logdir=$resultsdir/log

#iterate over all files in fastqdir, extract file identifier, run htseq-qa on file and write output into outputdir
for file in $fastqdir/*.fastq; do                                                                                                                    

    #if [[ $file =~ .*_2.fastq ]]  #make sure to exclude one  pair end read
    #then
    #    continue
    #fi

    name=`echo $file | perl -pe 's/.*\/(.*)\.fastq/$1/'`  
                                  
    if [ ! -f $resultsdir/${name}_htseq-qa_processed_rerun.pdf ]                                                                                                                   
    then                                                                                     
        echo $name
        qsub -N htseq-qa_processed_rerun_$name -l mem=$mem,walltime=3:00:00 -o $logdir/htseq-qa_processed_rerun_${name}.log -e $logdir/htseq-qa_processed_rerun_${name}.err -v outfile=$resultsdir/${name}_htseq-qa_processed_rerun\.pdf,file=$file htseq-qa_on_fastq.sh
    fi                                                                                                                                                 
done 
