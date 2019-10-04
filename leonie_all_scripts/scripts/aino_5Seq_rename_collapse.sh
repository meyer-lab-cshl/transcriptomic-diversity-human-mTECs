#!/bin/bash


samfile=$samfile
collapsed_file=$collapsed_file
mc_count_file=$mc_count_file
pos_count_file=$pos_count_file
collapse_info=$collapse_info

python /home/stroemic/hiwi_16/scripts/scripts_aino/5Seq_rename-and-collapse.py --samfile $samfile --out-collapsed $collapsed_file --out-molecule-count $mc_count_file --out-position-count $pos_count_file --out-info $collapse_info 
