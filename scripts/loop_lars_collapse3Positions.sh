#!/bin/bash

#example: bash loop_lars_collapse3Positions.sh /home/stroemic/hiwi_16/data/Summary_counts/mouse/ /home/stroemic/hiwi_16/analysis/lars_perl/mouse_consensus_sites.csv /home/stroemic/hiwi_16/analysis/lars_perl/mouse_consensus_sites.summary.csv 12 10  

#save inputdirectories into variables                                                                                                                                         
count_dir=$1
outfile_collapsed=$2
outfile_summarised=$3
max_distance=$4
min_no_bc=$5

mem=100gb    

out_dir=$(dirname "${outfile_collapsed}")
mkdir -p $out_dir/log/collapse_consensus

if [ -f $outfile_collapsed ]
then
    qsub -N collapse_consensus -l mem=$mem,walltime=36:00:00 -o $out_dir/log/collapse_consensus/collapse_consensus_mouse.log -e $out_dir/log/collapse_consensus/collapse_consensus_mouse.err -v count_dir=$count_dir,outfile_collapsed=$outfile_collapsed,outfile_summarised=$outfile_summarised,max_distance=$max_distance,min_no_bc=$min_no_bc  lars_collapse3Positions.sh
fi
