#!/bin/bash

file=$file
outdir=$outdir
mode=$mode
no_motifs=$no_motifs

/home/stroemic/software/miniconda2/envs/rnaseq/bin/meme $file -oc $outdir -dna -mod $mode -nmotifs $no_motifs 
