#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=2G
#$ -N SalmonTE
#$ -o SalmonTE_output.txt
#$ -e SalmonTE_output.txt

echo '########'
echo 'New run'
echo '########'

cat $1 | while read line || [[ -n $line ]];do

	SalmonTE.py quant \
	--reference=hs \
	--outpath=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/data/$line \
	--num_threads=4 \
	/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis/data/$line/fastp
	
	sleep 3h

done

  
  

