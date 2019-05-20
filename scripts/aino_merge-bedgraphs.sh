#!/bin/bash


xdir=$xdir
filePlus=$filePlus
fileMinus=$fileMinus
outfile=$outfile

echo 'before code'
echo $xdir
echo $filePlus
echo $fileMinus
echo $outfile

bash /home/stroemic/hiwi_16/scripts/scripts_aino/merge-bedgraphs.sh $xdir $filePlus $fileMinus $outfile
