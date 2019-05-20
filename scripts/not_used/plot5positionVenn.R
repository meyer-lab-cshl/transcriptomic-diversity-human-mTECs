# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      A simple script to plot Venn diagrams based on a in input table defining 5' positions.
#	(collapsed-counts-table_genomic-chr_single-nt.tsv)
#
# cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args
#	infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv'
#	outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-212_single-nt_all-counts.tiff'
#	samples='212_hi,212_lo'
#	cutoff=1
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

options(scipen=20)

library(VennDiagram)
library(RColorBrewer)

####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
	stop("No arguments provided!\n")
} else {

	infile =ifelse( length(grep('infile', args) )>0, strsplit(args[grep('infile', args)], '=')[[1]][2], stop('No infile argument provided!'))
	outfile =ifelse( length(grep('outfile', args) )>0, strsplit(args[grep('outfile', args)], '=')[[1]][2], stop('No outfile argument provided!'))
	samples =ifelse( length(grep('samples', args) )>0, strsplit(args[grep('samples', args)], '=')[[1]][2], stop('No samples argument provided!'))
	samples=strsplit(samples, ',')[[1]]
	cutoff =ifelse( length(grep('cutoff', args) )>0, strsplit(args[grep('cutoff=', args)], '=')[[1]][2], stop('No cutoff argument provided!'))
	cutoff=as.integer(cutoff)

}

####### PLOT

t = read.delim(infile)
t2=t[,unlist(lapply(samples, function(x) { grep(x, colnames(t)) }))]

positionLists = lapply( 1:length(t2[1,]), function(x) {
	return(rownames( t2[ which(t2[,x]>=cutoff), ] ))
})
names(positionLists)=sub('^X', '', colnames(t2))

# venn.diagram(list(ESC1=ESC1, ESC2=ESC2), outfile,
venn.diagram(positionLists, outfile,
	col=rep(NA, length(positionLists)),
	fill=brewer.pal(12, "Set3")[1:length(positionLists)], # c( rgb(0,.5,.1,.5), rgb(0,0,1,.5) ),
	force.unique=T,
	cex=1.5,
	label.col="black",
	cat.cex=2.5,
	cat.col=brewer.pal(12, "Set3")[1:length(positionLists)], # c( rgb(0,.5,.1,.5), rgb(0,0,1,.5) ),
	margin=0.1
)

