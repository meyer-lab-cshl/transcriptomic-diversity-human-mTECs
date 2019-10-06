#!/bin/bash

file=$file
background=$background
outdir=$outdir
no_motifs=$no_motifs

/home/stroemic/software/miniconda2/envs/rnaseq/bin/dreme -p $file -n $background -oc $outdir -dna -png #-m $no_motifs 
