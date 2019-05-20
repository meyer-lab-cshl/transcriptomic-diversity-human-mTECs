#!/bin/bash

bamfile=$bamfile
outfile=$outfile

/home/stroemic/software/miniconda2/envs/rnaseq/bin/genomeCoverageBed -5 -dz -ibam $bamfile > $outfile
