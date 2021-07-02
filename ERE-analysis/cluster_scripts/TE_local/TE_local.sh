#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=32G
#$ -N TE_local
#$ -o TE_local_output.txt
#$ -e TE_local_output.txt

cd $TE_HOME/data/RNA_seq/muscle/bam_files/TE_local_test_set                    

for FILE in *_Aligned.out.bam; do
  TElocal \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TElocal_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
  --project TE_local_${FILE}

done

cd $TE_HOME/data/RNA_seq/testis_jason/TE_local_test_set

for FILE in *_Aligned.out.bam; do
  TElocal \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TElocal_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
  --project TE_local_${FILE}

done

cd $TE_HOME/data/RNA_seq/ovaries_jason/TE_local_test_set

for FILE in *_Aligned.out.bam; do
  TElocal \
  -b ${FILE} \
  --GTF $TE_HOME/index/annotations/human.GRCh38.gtf \
  --TE $TE_HOME/index/annotations/TElocal_prebuilt_indices/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
  --project TE_local_${FILE}

done



