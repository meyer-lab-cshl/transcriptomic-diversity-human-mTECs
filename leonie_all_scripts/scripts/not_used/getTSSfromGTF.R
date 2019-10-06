# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      This processes a GTF file to write out an IRanges object of unique TSS positions
#
# cat /g/steinmetz/project/mTEC_Seq/src/getTSSfromGTF.R | R --slave --args
#	gtf='/g/steinmetz/genome/Homo_sapiens/37.68/annotation/gtf/Homo_sapiens.GRCh37.68.chrOnly.gtf'
#	outfile='/g/steinmetz/project/mTEC_Seq/annotation/TSS_IRanges.rda'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(IRanges)


####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
	stop("No arguments provided!\n")
} else {

	gtf =ifelse( length(grep('gtf', args) )>0, strsplit(args[grep('gtf=', args)], '=')[[1]][2], stop('No gtf argument provided!'))
	outfile =ifelse( length(grep('outfile', args) )>0, strsplit(args[grep('outfile=', args)], '=')[[1]][2], stop('No outfile argument provided!'))
}


######## FIXED PARAMETERS

genomic=c( '1','10','11','12','13','14','15','16','17','18','19','2','3','4','5','6','7','8','9','20','21', '22','X','Y' )


####### DEFINE TSSs

gtf=read.delim(gtf,header=F)
gtf_genomic = gtf[which( gtf$V1 %in% genomic ), ]

# Take exons
exons=gtf[which(gtf$V3=='exon'),]

# Extract transcript ids, add to table
txids = lapply( as.character(exons[,9]), function(x) strsplit( strsplit(x, '; ')[[1]] [ grep( 'transcript_id', strsplit(x, '; ')[[1]] ) ], ' ' )[[1]][2] )

t= as.data.frame(matrix(nrow=length(exons[,1]), ncol=0))
t$chr = as.character(exons[,1])
t$strand = as.character(exons[,7])
t$start = as.integer(as.character(exons[,4]))
t$end = as.integer(as.character(exons[,5]))
t$txid = unlist(txids)

# Separate by strand, extract start of TSS exon, and make a table
plus=t[which(t$strand=='+'),]
minus=t[which(t$strand=='-'),]
plus_TSS = by(plus$start, plus$txid, min)
minus_TSS = by(minus$end, minus$txid, max)

plus_t = as.data.frame(matrix(nrow=length(plus_TSS), ncol=0))
plus_t$txid = names(plus_TSS)
plus_t$TSS = as.integer(as.character(plus_TSS))
minus_t = as.data.frame(matrix(nrow=length(minus_TSS), ncol=0))
minus_t$txid = names(minus_TSS)
minus_t$TSS = as.integer(as.character(minus_TSS))

TSSs = rbind(plus_t, minus_t)

m = merge(unique( t[,c(1,2,5)] ), TSSs, by.x='txid', by.y='txid')
# [1] 91568     4

## Gene gene_id information
#geneids = unlist(lapply( as.character(exons[,9]), function(x) tail(strsplit( strsplit(x, '; ')[[1]] [ grep( 'gene_id', strsplit(x, '; ')[[1]] ) ], ' ' )[[1]], 1) ))
#gids = unlist(lapply(m$txid, function(x) { geneids[which(txids %in% x)] }))

# --> Make an IRanges
TSS_ir = RangedData(
	IRanges(
	start= m$TSS,
	end= m$TSS),
	space=paste( m$chr,'.', m$strand,sep=''),
	chr= m$chr,
	strand= m$strand,
	txid=m$txid #, geneid=gids
)

save(TSS_ir, file=outfile)

