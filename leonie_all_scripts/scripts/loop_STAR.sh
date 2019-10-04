#!/bin/bash

#example: bash loop_STAR.sh /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed /home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT /home/stroemic/genomes/human/85/DNA/STAR_genome_IVT
#save inputdirectories into variables                                                                                                                             
fastqdir=$1                                                                                                                                                
resultsdir=$2
genomedir=$3
mem=35gb
                                                                                                                                                          
mkdir -p $resultsdir                                                                                                                                      
mkdir -p $resultsdir/log                                                                                                                                      

logdir=$resultsdir/log

#iterate over all files in fastqdir, extract file identifier, run fatsqc on file and write output into outputdir, directly extract files from zip archiven
for file in $fastqdir/*.fastq; do                                                                                                                    

    if [[ $file =~ .*_2.fastq ]]  #make sure to exclude one  pair end read
    then
        continue
    fi

    #name=`echo $file | perl -pe 's/.*\/(.*)\_(1|2).fastq/$1/'`  
    name=`echo $file | perl -pe 's/.*\/(.*)\.fastq/$1/'`  
                                  
    if [ ! -f $resultsdir/${name}/${name}_Log.final.out ]                                                                                                                   
    then                                                                                                                                                
        echo $name
        mkdir -p $resultsdir/${name} 
        #qsub -N STAR_$name -l mem=$mem,walltime=24:00:00,nodes=1:ppn=4 -o $logdir/STAR_${name}.log -e $logdir/STAR_${name}.err -v thread_number=4,genomedir=$genomedir,read1=$fastqdir/${name}_1.fastq,read2=$fastqdir/${name}_2.fastq,outprefix=$resultsdir/$name/${name}_ STAR.sh
        qsub -N STAR_$name -l mem=$mem,walltime=24:00:00,nodes=1:ppn=4 -o $logdir/STAR_${name}.log -e $logdir/STAR_${name}.err -v thread_number=4,genomedir=$genomedir,read1=$file,outprefix=$resultsdir/$name/${name}_ STAR.sh
    fi                                                                                                                                                 
done 
