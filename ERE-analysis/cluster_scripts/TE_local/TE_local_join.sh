#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=5G
#$ -N TE_local_join
#$ -o TE_local_join_output.txt
#$ -e TE_local_join_output.txt

TE_HOME=/grid/meyer/home/mpeacey/thymus-epitope-mapping/ERE-analysis
cd $TE_HOME/data/

for TISSUE in *; do
	
	echo ${TISSUE}
	cd ${TISSUE}

	declare -a file_array=()

  	for TABLE in TE_local_*.cntTable; do
    		file_array+=($TABLE)
  	done

  	counter=0
  	while [ $((counter+1)) -lt ${#file_array[@]} ]; do

    		if [ $counter -eq 0 ]; then

      			join ${file_array[$counter]} ${file_array[$((counter+1))]} > TE_local_output_$((counter+1)) 

    		else

      			join TE_local_output_${counter} ${file_array[$((counter+1))]} > TE_local_output_$((counter+1))

    		fi 

    	counter=$((counter+1))

  	done

	mv TE_local_output_$((counter)) ../TE_local_${TISSUE}.cntTable
  	rm TE_local_output_*
  	cd ..

done



