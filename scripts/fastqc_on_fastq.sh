#!/bin/bash

file=$file
fastqcresultsdir=$outdir

fastqc --extract -o $fastqcresultsdir  $file
