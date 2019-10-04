#!/bin/bash

#example: qsub -N 5Seq_sort_and_trim -l mem=10gb,walltime=24:00:00 -o /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/log/C3V0VACXX_sort_trim_rev.log -e /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/log/C3V0VACXX_sort_trim_rev.err -v read1=/home/stroemic/hiwi_16/data/internal/5Seq/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz,read2=/home/stroemic/hiwi_16/data/internal/5Seq/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_2_sequence.txt.gz,designfile=/home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/5Seq_C3V0VACXX.design,infofile=/home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/5Seq_C3V0VACXX.info,outfile=/home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/5Seq_C3V0VACXX.out aino_5Seq_sort_trim.sh

forward=$read1
reverse=$read2
design=$designfile
info=$infofile
output=$outfile

python /home/stroemic/hiwi_16/scripts/scripts_aino/5seq_sort-and-trim_rev.py --forward $read1 --reverse $read2 --design $design --sample-barcode 6 --molecular-barcode 8 --trim-molecular-barcode --info $info > $output
