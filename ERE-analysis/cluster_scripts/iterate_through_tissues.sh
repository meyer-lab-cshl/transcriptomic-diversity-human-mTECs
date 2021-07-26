#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N iterate_through_tissues
#$ -o iterate_output.txt
#$ -e iterate_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/
cd $TE_HOME/data/

#for TISSUE in $1;do
#	
#	cd $TISSUE
#
#	for FILE in *.bam;do
#	
#		qsub /grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/cluster_scripts/bam_to_fastp.sh $FILE
#
#	done
#
#	sleep 0.5h
#
#	cd ..
#	
#done

cat $1 | while read line || [[ -n $line ]];do

	cd $line/
	
	for FILE in *.bam;do

		qsub /grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/cluster_scripts/bam_to_fastp.sh $FILE			

	done

	sleep 1h

	cd ..

done
