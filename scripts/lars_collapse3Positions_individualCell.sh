#!/bin/bash

sc_file=$sc_file
isoform_deffile=$isoform_deffile
out_file_isoforms=$out_file_isoforms
out_file_genes=$out_file_genes

echo $sc_file
echo $isoform_deffile
echo $out_file_isoforms
echo $out_file_genes

cat $sc_file | perl /home/stroemic/hiwi_16/scripts/scripts_lars/collapse3Positions_individualCell.pl  $isoform_deffile $out_file_isoforms $out_file_genes
