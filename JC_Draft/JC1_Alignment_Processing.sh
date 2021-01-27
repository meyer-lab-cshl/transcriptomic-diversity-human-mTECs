#!/bin/bash
#Starting with preprocessed files:
#        /grid/meyer/home/hmeyer/data/tss/human/5Pseq/deduplicated
# Convert sam to bam

project_directory=${HOME}/TSS
in_dir=/grid/meyer/home/hmeyer/data/tss/human/5Pseq/deduplicated/
out_dir=${project_directory}/5PSeq_Data/
subjects=(pt212-hi pt212-lo pt214-hi pt214-lo pt221-hi pt221-lo pt226-hi pt226-lo pt87-hi pt87-lo)
for subject in ${subjects[@]}
do
samtools view -b ${in_dir}${subject}_Aligned.sortedByCoord.dedup.unique.fwd.sam > ${out_dir}${subject}_Aligned.sortedByCoord.dedup.unique.fwd.bam
done

#Create necessary folders
dirs=(
5PSeq_Data
CTSS/Merged_CTSS_BED/temp
CTSS/CAGEr_out/DFs
CTSS/Paraclu_out/Paraclu_Cluster
CTSS/Paraclu_out/Paraclu_BED/temp
CTSS/Paraclu_out/Paraclu_Trim
CTSS/for_Paraclu
CTSS/HOMER/HOMER_Motifs/Motifs_Low
CTSS/HOMER/HOMER_Motifs/Motifs_High
CTSS/HOMER/HOMER_High_Low)
for dir in ${dirs[@]}
do
mkdir -p ${project_directory}/${dir}
done



# To arrive at prepocessed files:
# Adapted from Hannah's snakemake files

#!/bin/bash
#$ -t 1-10
# qsub Alignment_JC.sh
# Starting files from /grid/meyer/home/hmeyer/data/tss/human/5Pseq/fastq_filtered
#     cp *.paired.fastq ${data_directory}

#project_directory=${HOME}/TSS
#data_directory=${project_directory}/5PSeq_Data
#output_directory=${project_directory}/5PSeq_Aligned
#final_directory=${project_directory}/5PSeq_Processed_Data
#STAR_index=/grid/meyer/home/jacarter/Reference_Genome/STAR_Index/hg38

#Array job
#subjects=(pt212-hi pt212-lo pt214-hi pt214-lo pt221-hi pt221-lo pt226-hi pt226-lo pt87-hi pt87-lo)
#subject=${subjects[${SGE_TASK_ID}-1]}

#Alignment
#read1=${data_directory}/${subject}_1.paired.fastq
#read2=${data_directory}/${subject}_2.paired.fastq
#STAR --runThreadN 4 \
#    --runMode alignReads \
#    --genomeDir ${STAR_index} \
#    --readFilesIn ${read1} ${read2} \
#    --outReadsUnmapped Fastq \
#    --outSAMattrRGline ID:${subject} SM:${subject} \
#    --quantMode GeneCounts \
#    --outFileNamePrefix ${output_directory}/${subject}_

#samtools index ${output_directory}/${subject}_Aligned.sortedByCoord.out.bam

#umi_tools group \
#            -I ${output_directory}/${subject}_Aligned.sortedByCoord.out.bam \
#            --paired \
#            --group-out=${output_directory}/${subject}_Aligned.sortedByCoord.group.tsv \
#            --output-bam \
#            -S ${output_directory}/${subject}_Aligned.sortedByCoord.group.bam

#Dedup
#umi_tools dedup \
#            -I ${output_directory}/${subject}_Aligned.sortedByCoord.group.bam \
#            -L ${output_directory}/${subject}_dedup_extract.log \
#            -S ${output_directory}/${subject}_Aligned.sortedByCoord.dedup.bam \
#            --paired \
#            --output-stats=${output_directory}/${subject}_Aligned.sortedByCoord.stats_

#Keep only uniquely mapped reads
#samtools view -b -q 255 ${output_directory}/${subject}_Aligned.sortedByCoord.dedup.bam > ${output_directory}/${subject}_Aligned.sortedByCoord.dedup.unique.bam

#Keep only forwrard reads
#samtools view -h -f 0x40 ${output_directory}/${subject}_Aligned.sortedByCoord.dedup.unique.bam > ${final_directory}/${subject}_Aligned.sortedByCoord.dedup.unique.fwd.bam
