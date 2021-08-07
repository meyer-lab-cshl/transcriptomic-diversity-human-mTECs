#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N TE_local_wrapper
#$ -o TE_local_wrapper_output.txt
#$ -e TE_local_wrapper_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis
cd $TE_HOME/data/

for TISSUE in *; do
	
	echo ${TISSUE}
	cd ${TISSUE}

	for FILE in *Aligned.out.bam; do

		qsub /grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/cluster_scripts/TE_local/TE_local.sh ${TISSUE} ${FILE}
	
	done
	
	cd ..

	sleep 3h

done



