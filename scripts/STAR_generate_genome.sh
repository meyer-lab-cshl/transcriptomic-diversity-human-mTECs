#!/bin/bash


thread_number=$thread_number
genomedir=$genomedir
genomefastadir=$genomefastadir
gtf=$gtf
length=$length #should be max(Readlength)-1


#generate genome locally on computing node as a lot of temporary files are created
#also prevents long running time, cause no need to read files over network connecting 
cd $PBS_SCRATCH_DIR/$PBS_JOBID

STAR --runThreadN $thread_number --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $genomefastadir --sjdbGTFfile $gtf --sjdbOverhang $length
# copy files to the wanted folder
cp * $genomedir
