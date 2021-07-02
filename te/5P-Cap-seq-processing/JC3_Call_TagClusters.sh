#!/bin/bash
#$ -t 1-10
# qsub JC3_Call_TagClusters.sh

# Run Paraclu on each sample indepedently
in_dir="/grid/meyer/home/jacarter/TSS/CTSS/for_Paraclu/"
out_dir="/grid/meyer/home/jacarter/TSS/CTSS/Paraclu_out/"
cluster_dir=${out_dir}"Paraclu_Cluster/"
trim_dir=${out_dir}"Paraclu_Trim/"
bed_dir_temp=${out_dir}"Paraclu_BED/temp/"
bed_dir=${out_dir}"Paraclu_BED/"
paraclu_in="_TPM_for_Paraclu.txt"
paraclu_out="_TPM_paraclu.txt"
paraclu_trim="_TPM_paraclu_simplified.txt"
bed="_TPM_paraclu_simplified_bed.txt"
bed_open="_TPM_paraclu_simplified_bed_"

# Define minimum tag cluster expression=2TPM, maximum cluster length=20bp
min_TPM=2
max_length=20

#Array Job
subjects=(pt212_hi pt212_lo pt214_hi pt214_lo pt221_hi pt221_lo pt226_hi pt226_lo pt87_hi pt87_lo)
subject=${subjects[$SGE_TASK_ID-1]}

#Paraclu - tag clusters
paraclu ${min_TPM} ${in_dir}${subject}${paraclu_in} > ${cluster_dir}${subject}${paraclu_out}

#Filter Paraclu - apply max length, min expression thresholds
paraclu-cut.sh -l ${max_length} ${cluster_dir}${subject}${paraclu_out} > ${trim_dir}${subject}${paraclu_trim}

#Process paraclu output, export as BED format
awk '{print $1 "\t" $3 "\t" $4 "\t" $2}' ${trim_dir}${subject}${paraclu_trim} > ${bed_dir_temp}${subject}${bed_open}"1.txt"
sed -i "s/$/\t${subject}\t${subject}_/" ${bed_dir_temp}${subject}${bed_open}"1.txt"
awk '{ print $0, NR }' ${bed_dir_temp}${subject}${bed_open}"1.txt" > ${bed_dir_temp}${subject}${bed_open}"2.txt"
awk '{new_var=$6$7; print $0, new_var}' ${bed_dir_temp}${subject}${bed_open}"2.txt"  > ${bed_dir_temp}${subject}${bed_open}"3.txt"
awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $8 "\t" $4}' ${bed_dir_temp}${subject}${bed_open}"3.txt"  > ${bed_dir}${subject}${bed}
