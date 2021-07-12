#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N rename
#$ -o rename_output.txt
#$ -e rename_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/
cd $TE_HOME/data/

for TISSUE in *;do
  
  cd ${TISSUE}
  rm -r temp

  for FILE in *.bam;do
	
	echo ${FILE}
	PREFIX=`echo ${FILE} | cut -d'_' -f 1`
	PREFIX=`echo ${PREFIX} | cut -d'-' -f 2-5`
	MIDDLEFIX=`echo ${TISSUE} | tr _ -`
	mv ${FILE} ${PREFIX}_${MIDDLEFIX}_GTEX.Aligned.out.bam

  done
 
  cd ..
  
done

