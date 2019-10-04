#!/bin/bash


count_dir=$count_dir
outfile_collapsed=$outfile_collapsed
outfile_summarised=$outfile_summarised
max_distance=$max_distance
min_no_bc=$min_no_bc

#echo $count_dir
#echo $outfile_collapsed
#echo $outfile_summarised
#echo $max_distance
#echo $min_no_bc

perl /home/stroemic/hiwi_16/scripts/scripts_lars/collapse3Positions.pl  $count_dir $outfile_collapsed $outfile_summarised $max_distance $min_no_bc
