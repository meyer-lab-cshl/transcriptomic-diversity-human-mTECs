# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      Plotting accuracy of IVT 'TSS' identification.
#
# cat /g/steinmetz/project/mTEC_Seq/src/getIVTaccuracy.R | R --slave --args
#	infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv'
#	indir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
#	pattern='countsPerMol.tsv'
#	outPlot='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/5Seq-mapping-of-IVT-TSS.pdf'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library("reshape2", lib.loc="/g/steinmetz/jaerveli/bin/R/")
library("ggplot2", lib.loc="/g/steinmetz/jaerveli/bin/R/")



####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
	stop("No arguments provided!\n")
} else {
	infile =ifelse( length(grep('infile', args) )>0, strsplit(args[grep('infile=', args)], '=')[[1]][2], stop('No infile argument provided!'))	
	indir =ifelse( length(grep('indir', args) )>0, strsplit(args[grep('indir=', args)], '=')[[1]][2], stop('No indir argument provided!'))	
	pattern =ifelse( length(grep('pattern', args) )>0, strsplit(args[grep('pattern=', args)], '=')[[1]][2], stop('No pattern argument provided!'))	
	outPlot =ifelse( length(grep('outPlot', args) )>0, strsplit(args[grep('outPlot=', args)], '=')[[1]][2], stop('No outPlot argument provided!'))	
}


####### HELPER FUNCTIONS

## For each IVT, check what % of reads are within +/- 2 nt from TSS
getStartPercent = function(ivt) {
	p = round( 100 * ( length( which( countsPerMol[which(countsPerMol$chr==ivt),]$pos >= 15
		& countsPerMol[which(countsPerMol$chr==ivt),]$pos <= 19 ) ) ) / ( length( countsPerMol[which(countsPerMol$chr==ivt),1] ) ), digits = 2 )
	return(p)
}


####### READ INPUT

files = list.files(path=indir, pattern=pattern)

t=as.data.frame(matrix(nrow=3, ncol=length(files)))
colnames(t)=unlist(lapply(files, function(x) strsplit(x, '_')[[1]][1] ))
rownames(t)=c('lys','phe','thr')

counter=1
for( f in files ) {
	countsPerMol = read.delim(paste(indir,f,sep=''))
	t[1:3,counter] = c( getStartPercent('pGIBS-LYS'), getStartPercent('pGIBS-PHE'), getStartPercent('pGIBS-THR') )
	counter=counter+1
}

### PLOT

mt = as.data.frame( matrix( nrow=length(rownames(t))*length(colnames(t)), ncol=0 ) )
mt$sample=rep( rownames(t), length(colnames(t))) # this is each bar.
mt$group=rep(colnames(t), each=3)
mt$ivt=rep( c('lys','phe', 'the'), length(colnames(t)))
mt$percent=unlist( lapply(t, c)) # c( t[,1], t[,2], t[,3], t[,4] ) # all 

pdf(outPlot)
ggplot() +
  geom_bar(data=mt, aes(y = percent, x = sample, fill = ivt), stat="identity",
           position='stack') +
  theme_bw() + 
  facet_grid( ~ group) + labs(title="Mapped 5 position +/- 2bp from known IVT TSS")

dev.off()

