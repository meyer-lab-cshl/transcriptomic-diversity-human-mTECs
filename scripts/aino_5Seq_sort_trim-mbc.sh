#!/bin/bash

#example: qsub -N 5Seq_sort_and_trim-mbc -l mem=2gb,walltime=6:00:00 -o /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/log/C83FLACXX_pt221-lo_sort_trim.log -e /home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/log/C83FLACXX_pt221-lo_sort_trim.err -v read1=/home/stroemic/hiwi_16/data/internal/5Seq/C83FLACXX_pt221-lo_16s000847-1-1_Clauder-Muenster_lane816s000847_1_sequence.txt.gz,read2=/home/stroemic/hiwi_16/data/internal/5Seq/C83FLACXX_pt221-lo_16s000847-1-1_Clauder-Muenster_lane816s000847_2_sequence.txt.gz,designfile=/home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/5Seq_C83FLACXX_pt221-lo.design,infofile=/home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/5Seq_C83FLACXX_pt221-lo.info,outfile=/home/stroemic/hiwi_16/data/internal/5Seq/5Seq_fastq_processed/5Seq_C83FLACXX_pt221-lo.out aino_5Seq_sort_trim-mbc.sh

forward=$read1
reverse=$read2
design=$designfile
info=$infofile
output=$outfile

python /home/stroemic/hiwi_16/scripts/scripts_aino/5seq_sort-and-trim-mbc.py --forward $read1 --reverse $read2 --design $design --sample-barcode 6 --molecular-barcode 8 --trim-molecular-barcode --info $info > $output
