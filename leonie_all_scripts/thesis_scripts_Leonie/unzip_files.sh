#!/bin/bash


#example: bash unzip_files.sh /home/stroemic/Data/RNASeq
dir=$1
mem=2000
mkdir -p /home/stroemic/Data/log

logdir=/home/stroemic/Data/log
#PBS -l walltime=1:00:00
#PBS -l mem=1024m
#PBS -M stroemic@dkfz.de
#PBS -o /home/stroemic/Data/log/unzip/unzip.log

for file in $dir/*.fastq.gz; do
    name=`echo $file | perl -pe 's/.*\/(.*)\.fastq\.gz/$1/'` 
    echo $name
    if [ ! -f $dir/$name.fastq ]
    then
        qsub -N unzip_$name -l mem=$mem,walltime=1:00:00 -o $logdir/unzip/unzip_${name}.log -e $logdir/unzip/unzip_${name}.err "java -jar ~/software/CRAM/cramtools-3.0.jar fastq -I $file --reverse > $fastqdir/$name.fastq"
    fi
done 
