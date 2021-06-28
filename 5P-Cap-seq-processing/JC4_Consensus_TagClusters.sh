#!/bin/bash

in_dir="/grid/meyer/home/jacarter/TSS/CTSS/Paraclu_out/Paraclu_BED/"
tmp_dir="/grid/meyer/home/jacarter/TSS/CTSS/Merged_CTSS_BED/temp/"
out_dir="/grid/meyer/home/jacarter/TSS/CTSS/Merged_CTSS_BED/"
bed="_TPM_paraclu_simplified_bed.txt"
union_bed="All_Thymus_bed.bed"
sorted_union_bed="All_Thymus_sorted_bed.bed"
merged_bed="All_thymus_merged.txt"

bedops --everything \
      ${in_dir}{"pt212_hi","pt212_lo","pt214_hi","pt214_lo","pt221_hi","pt221_lo","pt226_hi","pt226_lo","pt87_hi","pt87_lo"}${bed} \
      > ${tmp_dir}${union_bed}

sort-bed ${tmp_dir}${union_bed} > ${tmp_dir}${sorted_union_bed}

bedtools merge -i ${tmp_dir}${sorted_union_bed} -s -d 20 -c 4,5,6 \
      -o distinct,distinct,distinct > ${out_dir}${merged_bed}
