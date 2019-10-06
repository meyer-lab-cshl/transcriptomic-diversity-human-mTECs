# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 28 March 2013
# Description:
#      Preparing a gene annotation objects for assignment.
#       All exons, plus a 200bp extension upstream of TSS.
#       Note that in the output IR object, the extended regions have a biotype 'extension'.
#               Thus, they can be easily de-selected from the data if one does not want to use extension
#               (or wants an extension of a different length).
#
# cat /g/steinmetz/project/mTEC_Seq/src/prepareGeneAssignmentObject.R | R --slave --args
#       gtf='/g/steinmetz/genome/Homo_sapiens/37.68/annotation/gtf/Homo_sapiens.GRCh37.68.chrOnly.gtf'
#       outfile='/g/steinmetz/project/mTEC_Seq/annotation/genes-by-exons-plus-200-TSS-extension_IRanges.rda'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
        stop("No arguments provided!\n")
} else {
        gtf =ifelse( length(grep('gtf', args) )>0, strsplit(args[grep('gtf=', args)], '=')[[1]][2], stop('No gtf argument provided!'))
        outfile =ifelse( length(grep('outfile', args) )>0, strsplit(args[grep('outfile=', args)], '=')[[1]][2], stop('No outfile argument provided!'))
}


####### FIXED PARAMETERS

genomic=c( '1','10','11','12','13','14','15','16','17','18','19','2','3','4','5','6','7','8','9','20','21','22','X','Y' )

library(IRanges)


#### MAKE A GENE OBJECT (all exons of a gene + 200bp added to the annotated TSS)

anno = read.delim(gtf, header=F)

exons = anno[which(anno$V3=='exon'),]

geneids = unlist(lapply( as.character(exons$V9), function(x) { strsplit( strsplit(x, '; ')[[1]][grep( 'gene_id', strsplit(x, '; ')[[1]] )], 'gene_id ' )[[1]][2] }))
txids = unlist(lapply( as.character(exons$V9), function(x) { strsplit( strsplit(x, '; ')[[1]][grep( 'transcript_id', strsplit(x, '; ')[[1]] )], 'transcript_id ' )[[1]][2] }))
biotype = unlist(lapply( as.character(exons$V9), function(x) { strsplit( strsplit(x, '; ')[[1]][grep( 'gene_biotype', strsplit(x, '; ')[[1]] )], 'gene_biotype ' )[[1]][2] }))

## Collapse so that the exact same exons do not occur multiple times.
t = cbind(exons[,c(1,4,5,7)], geneids, txids, biotype) # 
t2 = unique(t[,c(1,2,3,4,5,7)])


### Extend each TX by 200bp upstream of its TSS

plus= t[which(t$V7=='+'),]
minus= t[which(t$V7=='-'),]

plus_start = by(plus$V4, as.character(plus$txids), min)
minus_start = by(minus$V5, as.character(minus$txids), max)

t3=as.data.frame(matrix(nrow=length(t$geneids), ncol=0))
t3$txid=as.character(t$txids)
t3$geneid=as.character(t$geneids)
t3$chr=as.character(t$V1)
t3=unique(t3)

x=as.data.frame(matrix(nrow=length(names(plus_start)), ncol=0))
x$txid=names(plus_start)
x$start = as.integer(plus_start) - 200
x$end = as.integer(plus_start)
x$strand = rep('+', length(x$txid ))
x$biotype = rep('extension', length(x$txid ))
plus_extra = merge(x, t3, by.x='txid', by.y='txid')
plus_extra = unique(plus_extra[, c(2,3,4,5,6,7)])

x=as.data.frame(matrix(nrow=length(names(minus_start)), ncol=0))
x$txid=names(minus_start)
x$start = as.integer(minus_start)
x$end = as.integer(minus_start) + 200
x$strand = rep('-', length(x$txid ))
x$biotype = rep('extension', length(x$txid ))
minus_extra = merge(x, t3, by.x='txid', by.y='txid')
minus_extra = unique(minus_extra[, c(2,3,4,5,6,7)])

### Make an IRanges object.

gene_IR = RangedData(
        IRanges(
        start= c( t2$V4, plus_extra$start, minus_extra$start),
        end= c( t2$V5, plus_extra$end, minus_extra$end )),
        space= c( paste( t2$V1,'.', t2$V7,sep=''), paste( plus_extra$chr, '.', plus_extra$strand,sep=''), paste( minus_extra$chr, '.', minus_extra$strand,sep='')),
        chr= c( as.character(t2$V1), as.character(plus_extra$chr), as.character(minus_extra$chr) ),
        strand= c( as.character(t2$V7), as.character(plus_extra$strand), as.character(minus_extra$strand) ),
        id= c( as.character(t2$geneids), as.character(plus_extra$geneid), as.character(minus_extra$geneid) ),
        biotype= c( as.character(t2$biotype), as.character(plus_extra$biotype), as.character(minus_extra$biotype) )
)

gene_IR = gene_IR[which(gene_IR$chr %in% genomic),]

save(gene_IR, file=outfile)

