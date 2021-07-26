#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N SalmonTE
#$ -o SalmonTE_output.txt
#$ -e SalmonTE_output.txt

echo '########'
echo 'New run'
echo '########'

echo $2
cat $1 | while read line || [[ -n $line ]];do
	
#	if [ $2=='raw' ]
#	then
#		
#		echo 'Calculating raw counts for tissue' $line
#
#		SalmonTE.py quant \
 #       	--reference=hs \
  #      	--outpath=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/data/$line/raw_Salmon_counts/ \
   #     	--exprtype=count \
#		--num_threads=4 \
 #       	/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/data/$line/fastp
#	fi
	if [ $2=='TPM' ]
	then
		
		echo 'Calculating normalized counts for tissue' $line

		SalmonTE.py quant \
		--reference=hs \
		--outpath=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/data/$line \
		--num_threads=4 \
		/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/data/$line/fastp
	
	fi

	sleep 20m

done
