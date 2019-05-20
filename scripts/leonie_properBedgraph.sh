#!/bin/bash

xdir=$xdir
name=$name
bedgraph=$bedgraph
outfile=$outfile

python /home/stroemic/hiwi_16/scripts/properBedgraph.py $bedgraph $outfile

echo 'finished python'
echo $xdir
echo $name
echo $bedgraph
echo $outfile

echo 'track type=bedGraph name="" description="BedGraph format" visibility=full color=200,100,0' > $xdir/descript.txt
cat $xdir/descript.txt $outfile > $xdir/${name}_5p.dz-collapsed.upload.bedgraph
