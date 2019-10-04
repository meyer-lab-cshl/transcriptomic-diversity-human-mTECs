# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 28 March 2013
# Description:
#      Assignment of 5' data to annotation.
#	This script is (in principle) a general assignment script, but it was made and tested for gene annotation.
#
# cat /g/steinmetz/project/mTEC_Seq/src/assing2genes.R | R --slave --args
#	data='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv'
#	assignment='/g/steinmetz/project/mTEC_Seq/annotation/genes-by-exons-plus-200-TSS-extension_IRanges.rda'
#	samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo'
#	outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/5Seq-at-genes.tsv'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(IRanges)

####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
        stop("No arguments provided!\n")
} else {
        data =ifelse( length(grep('data', args) )>0, strsplit(args[grep('data=', args)], '=')[[1]][2], stop('No data argument provided!'))
        assignment =ifelse( length(grep('assignment', args) )>0, strsplit(args[grep('assignment=', args)], '=')[[1]][2], stop('No assignment argument provided!'))
        samples =ifelse( length(grep('samples', args) )>0, strsplit(args[grep('samples=', args)], '=')[[1]][2], stop('No samples argument provided!'))
        samples=strsplit(samples, ',')[[1]]
        outfile =ifelse( length(grep('outfile', args) )>0, strsplit(args[grep('outfile=', args)], '=')[[1]][2], stop('No outfile argument provided!'))
}


####### HARDCODED PARAMETERS

genomic=c( '1','10','11','12','13','14','15','16','17','18','19','2','3','4','5','6','7','8','9','20','21','22','X','Y' )


####### IR OUT OF 5SEQ DATA

fiveSeq = read.delim(data)
fiveSeq_genomic = fiveSeq[which(fiveSeq$chr %in% genomic), ]

fiveSeq_genomic_ir = RangedData(
        IRanges(
        start= fiveSeq_genomic$pos,
        end= fiveSeq_genomic$pos),
        space=paste( fiveSeq_genomic$chr,'.', fiveSeq_genomic$strand,sep=''),
        chr= fiveSeq_genomic$chr,
        strand= fiveSeq_genomic$strand,
        id=paste(fiveSeq_genomic$chr, fiveSeq_genomic$strand, fiveSeq_genomic$pos,sep='_'),
        fiveSeq_genomic[, 4:length(fiveSeq_genomic[1,])]
)


####### ASSIGN

load(assignment)

fiveSeq_at_genes = findOverlaps(fiveSeq_genomic_ir, gene_IR)
mfo = as.matrix(fiveSeq_at_genes)
oId=mfo[,1]
tId=mfo[,2]
gene_assignment =data.frame(
	 dataID=oId,
	 annoID=tId,
        gene_id= gene_IR$id[tId],
        #chr= gene_IR$chr[tId],
        #strand= gene_IR$strand[tId],
        fiveSeq_id = fiveSeq_genomic_ir$id[oId],
	 as.data.frame(fiveSeq_genomic_ir[oId, unlist( lapply(samples, function(x) { grep(x, colnames(fiveSeq_genomic_ir)) })) ])[,5:(4+length(samples))]
)

gene_assignment2 = unique( gene_assignment[,3:(4+length(samples))] )
countsPerGene = by(gene_assignment2[,3:(2+length(samples))], as.character(gene_assignment2$gene_id), colSums) # Approx 1 minute.
t = do.call(rbind, countsPerGene)

write.table(t, file=outfile, quote=F, sep='\t')

