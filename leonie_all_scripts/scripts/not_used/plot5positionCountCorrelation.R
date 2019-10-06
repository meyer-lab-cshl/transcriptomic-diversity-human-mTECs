# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      A simple script to plot 5' position count correlation for used specified samples and input table
#	(collapsed-counts-table_genomic-chr_single-nt.tsv)
#
# cat /g/steinmetz/project/mTEC_Seq/src/plot5positionCountCorrelation.R | R --slave --args
#	infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv'
#	outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/5pos-count-correlation.tiff'
#	samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

options(scipen=20)

####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
	stop("No arguments provided!\n")
} else {

	infile =ifelse( length(grep('infile', args) )>0, strsplit(args[grep('infile', args)], '=')[[1]][2], stop('No infile argument provided!'))
	outfile =ifelse( length(grep('outfile', args) )>0, strsplit(args[grep('outfile', args)], '=')[[1]][2], stop('No outfile argument provided!'))
	samples =ifelse( length(grep('samples', args) )>0, strsplit(args[grep('samples', args)], '=')[[1]][2], stop('No samples argument provided!'))
	samples=strsplit(samples, ',')[[1]]

}

t=read.delim(infile)
t2=t[,unlist(lapply(samples, function(x) { grep(x, colnames(t)) }))]
colnames(t2)=sub('_', ' ', sub('^X', '', colnames(t2)) )

#pdf(outfile)
tiff(outfile, width = 600, height = 600)
pairs(log2(t2), pch=16, cex=0.6, col="#00000010", main="5' position count correlation (single nt) (log2)")
dev.off()

