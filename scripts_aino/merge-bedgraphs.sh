#!/bin/bash
 
# A shell script to merge the output of plus and minus single bp bedgraphs
#dir="/g/steinmetz/project/PROMPT-seq/run/2013-05-17_C235KACXX_Ntini_pelechan03_13s003843/collapsed"
#plus="plus.bg"
#minus="minus.bg"
#output="output.bg"
 
#Leonie: Input already full path; deleted $dir before files

xdir=$1
plus=$2
minus=$3
output=$4
awk '{print $0, "+"}' $plus > $xdir/tx_plus.bg
awk '{print $0, "-"}' $minus > $xdir/tx_minus.bg
cat $xdir/tx_plus.bg $xdir/tx_minus.bg > $xdir/tx.bg
sort -k1,1 -k2,2n < $xdir/tx.bg > $xdir/tx3.bg
awk '{OFS="\t"; print $1,$2,$3,NR,$4,$5}' $xdir/tx3.bg > $output
rm $xdir/tx_plus.bg
rm $xdir/tx_minus.bg
rm $xdir/tx.bg
rm $xdir/tx3.bg
