# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 28 March 2013
# Description:
#	
#	This script identifies differential expression using DESeq.
#	Note: the order of samples (as defined by samples argument) and group ids (groupIDs) must be the correct one:
#	The script labels the columns accordingly.
#
# cat /g/steinmetz/project/mTEC_Seq/src/getDifferentiallyExpressedGenes.R | R --slave --args
#	geneTable='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/5Seq-at-genes.tsv'
#	outSig='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low-significant.tsv'
#	outAll='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low-all.tsv'
#	outPlot='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low-MA.png'
#	outInfo='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low.info'
#	sigCO=0.05
#	samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo'
#	groupIDs='hi,lo,hi,lo,hi,lo'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

# print(args)

if(length(args)==0){
        stop("No arguments provided!\n")
} else {
        geneTable =ifelse( length(grep('geneTable', args) )>0, strsplit(args[grep('geneTable=', args)], '=')[[1]][2], NA)
	 if(is.na(geneTable)) { stop('No geneTable argument provided!') } 
        outSig =ifelse( length(grep('outSig', args) )>0, strsplit(args[grep('outSig=', args)], '=')[[1]][2], NA)
	 if(is.na(outSig)) { stop('No outSig argument provided!') } 
        outAll =ifelse( length(grep('outAll', args) )>0, strsplit(args[grep('outAll=', args)], '=')[[1]][2], NA)
	 if(is.na(outAll)) { stop('No outAll argument provided!') } 
        outPlot =ifelse( length(grep('outPlot', args) )>0, strsplit(args[grep('outPlot=', args)], '=')[[1]][2], NA)
	 if(is.na(outPlot)) { stop('No outPlot argument provided!') } 
        outInfo =ifelse( length(grep('outPlot', args) )>0, strsplit(args[grep('outPlot=', args)], '=')[[1]][2], NA)
	 if(is.na(outInfo)) { stop('No outInfo argument provided!') } 
        samples =ifelse( length(grep('samples', args) )>0, strsplit(args[grep('samples=', args)], '=')[[1]][2], NA)
	 if(is.na(samples)) { stop('No samples argument provided!') } 
        samples=strsplit(samples, ',')[[1]]
        groupIDs =ifelse( length(grep('groupIDs', args) )>0, strsplit(args[grep('groupIDs=', args)], '=')[[1]][2], NA)
	 if(is.na(groupIDs)) { stop('No groupIDs argument provided!') } 
        groupIDs =strsplit(groupIDs, ',')[[1]]
        sigCO =ifelse( length(grep('sigCO', args) )>0, strsplit(args[grep('sigCO=', args)], '=')[[1]][2], NA)
	 if(is.na(sigCO)) { stop('No sigCO argument provided!') } 
	 sigCO=as.numeric(sigCO)
}


library(DESeq2)
library(biomaRt)


####### RUN DESEQ2 ON THE DATA

t = read.delim(geneTable)

# Order the table to be in input order.
t = t[ , unlist( lapply(samples, function(x) { grep(x, colnames(t)) }) )]

myColData=as.data.frame(matrix(nrow=length(samples), ncol=0))
myColData$condition = as.factor(groupIDs) # as.factor(c('high','low','high','low','high','low'))

rownames(myColData)=colnames(t)

dds <- DESeqDataSetFromMatrix(countData = t, colData = myColData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

## Write out table of differential expression results.
write.table(as.data.frame(res), file=outAll, quote=F, sep='\t', row.names=F)


### ATTACH INFORMATION TO THE TABLE 
# gene name, description, and coordinates (using bioMART).

sig = as.data.frame(res[which(res$padj<sigCO),])

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene2mart = getBM(attributes = c("ensembl_gene_id", "external_gene_id", "description", "chromosome_name", "strand", "start_position", "end_position"),
		filters = c("ensembl_gene_id"),
		values = rownames(sig), mart = mart)

gene2mart$description2 = unlist(lapply(gene2mart$description, function(x) { strsplit( x, " [", fixed=T)[[1]][1] }))
rownames(gene2mart)=gene2mart$ensembl_gene_id

y = merge( sig[,c(2,6)], gene2mart[,c(2,8,4:7)], by.x="row.names", by.y="row.names" ) #'ensembl_gene_id' )
y2=y[,c(1,4,2,3,5, 6:9)]
t = y2[order(y2$padj),]
colnames(t)[c(1,2,5)]=c('gene_id', 'gene_name', 'description')

### WRITE OUT THE TABLE
write.table(t, file=outSig, quote=F, sep='\t', row.names=F)


### PLOT
### Make the DESeq2 MA plot.
png(outPlot)
plotMA(dds)
dev.off()

