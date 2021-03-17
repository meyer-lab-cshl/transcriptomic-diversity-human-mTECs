#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=0.5G
#$ -N TE_local
#$ -o TE_local_output.txt
#$ -e TE_local_output.txt

cd $TE_HOME/data

for FILE in *.fastq_Aligned.out.bam; do

  cat <<EOF TE_local_${FILE}.sh
  #!/bin/bash

  EOF
  
done
