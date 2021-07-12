#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N copy
#$ -o copy_output.txt
#$ -e copy_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/
cd $TE_HOME/data/

cd $GTEX

for TISSUE in *;do

	cp -r ${TISSUE} $TE_HOME/data

done
