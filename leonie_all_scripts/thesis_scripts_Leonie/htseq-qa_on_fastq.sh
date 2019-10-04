#!/bin/bash

file=$file
outfile=$outfile

htseq-qa --type=fastq --outfile=$outfile $file

