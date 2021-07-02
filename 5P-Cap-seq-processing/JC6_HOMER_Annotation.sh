#!/bin/bash

export PATH="$HOME/ToolKit/HOMER/bin:$PATH"

in_dir="/grid/meyer/home/jacarter/TSS/CTSS/HOMER/"
out_dir="/grid/meyer/home/jacarter/TSS/CTSS/HOMER/HOMER_Out/"
annotatePeaks.pl ${in_dir}forHomer.txt hg38 -size -given > ${out_dir}Homerout.txt
annotatePeaks.pl ${in_dir}forHomer.txt hg38 -CpG -size -500,100  > ${out_dir}Homerout_CpG.txt
annotatePeaks.pl ${in_dir}forHomer.txt hg38 -size -2500,2500 -hist 150 -di > ${out_dir}Homerout_Nuc.txt
annotatePeaks.pl ${in_dir}forHomer.txt hg38 -size -200,200 -hist 1 -m ${HOME}/ToolKit/HOMER/motifs/tata.motif >  ${out_dir}Homer_TATA.txt

out_dir2="/grid/meyer/home/jacarter/TSS/CTSS/HOMER/HOMER_High_Low/"
annotatePeaks.pl ${in_dir}Homer_High.txt hg38 -size -2500,2500 -hist 150 -di > ${out_dir2}Homerout_Nuc_High.txt
annotatePeaks.pl ${in_dir}Homer_Low.txt hg38 -size -2500,2500 -hist 150 -di > ${out_dir2}Homerout_Nuc_Low.txt
annotatePeaks.pl ${in_dir}Homer_High.txt hg38 -size -200,200 -hist 1 -m ${HOME}/ToolKit/HOMER/motifs/tata.motif >  ${out_dir2}Homer_TATA_High.txt
annotatePeaks.pl ${in_dir}Homer_Low.txt hg38 -size -200,200 -hist 1 -m ${HOME}/ToolKit/HOMER/motifs/tata.motif >  ${out_dir2}Homer_TATA_Low.txt

out_dir3="/grid/meyer/home/jacarter/TSS/CTSS/HOMER/HOMER_Motifs/"
findMotifsGenome.pl  ${in_dir}Homer_High.txt hg38 ${out_dir3}Motifs_High
findMotifsGenome.pl  ${in_dir}Homer_Low.txt hg38 ${out_dir3}Motifs_Low
