# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      Preparing various annotation objects for assignment.
#
# cat /g/steinmetz/project/mTEC_Seq/src/makeAssignmentPies.R | R --slave --args
#	data='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv'
#	assignment='/g/steinmetz/project/mTEC_Seq/annotation/assignment_IRanges.rda'
#	extension=100
#	samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo'
#	outPdfPosition='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/plots/5Seq-assignment-positions_TSS-100.pdf'
#	outPdfCount='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/plots/5Seq-assignment-counts_TSS-100.pdf'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


library(IRanges)
library(RColorBrewer)

####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
	stop("No arguments provided!\n")
} else {
	data =ifelse( length(grep('data', args) )>0, strsplit(args[grep('data=', args)], '=')[[1]][2], stop('No data argument provided!'))
	assignment =ifelse( length(grep('assignment', args) )>0, strsplit(args[grep('assignment=', args)], '=')[[1]][2], stop('No assignment argument provided!'))
	extension=ifelse( length(grep('extension', args) )>0, strsplit(args[grep('extension=', args)], '=')[[1]][2], stop('No extension argument provided!'))
	extension=as.numeric(extension)
	samples =ifelse( length(grep('samples', args) )>0, strsplit(args[grep('samples=', args)], '=')[[1]][2], stop('No samples argument provided!'))
	samples=strsplit(samples, ',')[[1]]
	outPdfPosition=ifelse( length(grep('outPdfPosition', args) )>0, strsplit(args[grep('outPdfPosition=', args)], '=')[[1]][2], stop('No outPdfPosition argument provided!'))
	outPdfCount=ifelse( length(grep('outPdfCount', args) )>0, strsplit(args[grep('outPdfCount=', args)], '=')[[1]][2], stop('No outPdfCount argument provided!'))
}


####### FIXED PARAMETERS

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


####### LOAD AND EXTEND RELEVANT ANNOTATION

# loads: TSS_ir, exon_ir, gene_ir, gene_antisense_ir, fiveUTR_ir, threeUTR_ir, biprom_ir, postTTS_ir
load(assignment)

TSS_ir2= TSS_ir
start(TSS_ir2)=start(TSS_ir2)-extension
end(TSS_ir2)=end(TSS_ir2)+extension
TSS_ir= TSS_ir2

####### ASSIGN


## TSS

fiveSeq_within_TSS = findOverlaps(fiveSeq_genomic_ir, TSS_ir)
mfo = as.matrix(fiveSeq_within_TSS)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_TSS=data.frame(query=oId,subject=tId,
	TSS_id= TSS_ir$id[tId],
	fiveSeq_id = fiveSeq_genomic_ir$id[oId]
)
## What are the ones that DO have overlaps??
fiveSeq_within_TSS = fiveSeq_genomic_ir[ which( fiveSeq_genomic_ir$id %in% minfo_withOL_TSS$fiveSeq_id ), ]
## What are the ones that do NOT have overlaps?
fiveSeq_not_within_TSS = fiveSeq_genomic_ir[ which( !( fiveSeq_genomic_ir$id %in% minfo_withOL_TSS$fiveSeq_id) ), ]


## IF NOT TSS, 5UTR (fiveUTR_ir)

oltif_within_5UTR = findOverlaps(fiveSeq_not_within_TSS, fiveUTR_ir)
mfo = as.matrix(oltif_within_5UTR)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_5UTR=data.frame(query=oId,subject=tId,
	tx_id= fiveUTR_ir$txid[tId],
	fiveSeq_id = fiveSeq_not_within_TSS $id[oId]
)
## What are the ones that DO have overlaps??
fiveSeq_within_5UTR = fiveSeq_not_within_TSS[ which( fiveSeq_not_within_TSS$id %in% minfo_withOL_5UTR$fiveSeq_id ), ]
## What are the ones that do NOT have overlaps?
fiveSeq_not_within_5UTR = fiveSeq_not_within_TSS[ which( !( fiveSeq_not_within_TSS$id %in% minfo_withOL_5UTR$fiveSeq_id) ), ]


## TSS > other 5UTR > 3UTR > exon > intron > bi-promoter > antisense > other (intergenic)

oltif_within_3UTR = findOverlaps(fiveSeq_not_within_5UTR, threeUTR_ir)
mfo = as.matrix(oltif_within_3UTR)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_3UTR=data.frame(query=oId,subject=tId,
	tx_id= threeUTR_ir$txid[tId],
	fiveSeq_id = fiveSeq_not_within_5UTR$id[oId]
)
## What are the ones that DO have overlaps??
fiveSeq_within_3UTR = fiveSeq_not_within_5UTR[ which( fiveSeq_not_within_5UTR$id %in% minfo_withOL_3UTR$fiveSeq_id ), ]
## What are the ones that do NOT have overlaps?
fiveSeq_not_within_3UTR = fiveSeq_not_within_5UTR[ which( !( fiveSeq_not_within_5UTR$id %in% minfo_withOL_3UTR$fiveSeq_id) ), ]



## If not TSS, exon --> IF NOT 3UTR, EXONIC

oltif_within_exon = findOverlaps(fiveSeq_not_within_3UTR, exon_ir)
mfo = as.matrix(oltif_within_exon)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_exon=data.frame(query=oId,subject=tId,
	#exon_id= exon_ir$id[tId],
	fiveSeq_id = fiveSeq_not_within_3UTR$id[oId]
)
## What are the ones that DO have overlaps??
fiveSeq_within_exon = fiveSeq_not_within_3UTR[ which( fiveSeq_not_within_3UTR$id %in% minfo_withOL_exon$fiveSeq_id ), ]
## What are the ones that do NOT have overlaps?
fiveSeq_not_within_exon = fiveSeq_not_within_3UTR[ which( !( fiveSeq_not_within_3UTR$id %in% minfo_withOL_exon$fiveSeq_id) ), ]



## If no exon, intron (gene)

oltif_within_gene = findOverlaps(fiveSeq_not_within_exon, gene_ir)
mfo = as.matrix(oltif_within_gene)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_gene=data.frame(query=oId,subject=tId,
	gene_id = gene_ir$gene[tId],
	fiveSeq_id= fiveSeq_not_within_exon$id[oId]
)

fiveSeq_within_intron = fiveSeq_not_within_exon[ which( fiveSeq_not_within_exon$id %in% minfo_withOL_gene$fiveSeq_id ), ]
fiveSeq_not_within_intron = fiveSeq_not_within_exon[ which( !( fiveSeq_not_within_exon$id %in% minfo_withOL_gene$fiveSeq_id) ), ]


### TSS > other 5UTR > 3UTR > exon > intron > bi-promoter > antisense > after TTS > other (intergenic)

oltif_within_biprom = findOverlaps(fiveSeq_not_within_intron, biprom_ir)
mfo = as.matrix(oltif_within_biprom)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_biprom=data.frame(query=oId,subject=tId,
	id = biprom_ir$id[tId],
	fiveSeq_id= fiveSeq_not_within_intron$id[oId]
)

fiveSeq_within_biprom = fiveSeq_not_within_intron[ which( fiveSeq_not_within_intron$id %in% minfo_withOL_biprom$fiveSeq_id ), ]
fiveSeq_not_within_biprom = fiveSeq_not_within_intron[ which( !( fiveSeq_not_within_intron$id %in% minfo_withOL_biprom$fiveSeq_id) ), ]



### TSS > other 5UTR > 3UTR > exon > intron > bi-promoter > after TTS > antisense > other (intergenic)

oltif_within_biprom = findOverlaps(fiveSeq_not_within_biprom, postTTS_ir)
mfo = as.matrix(oltif_within_biprom)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_postTTS=data.frame(query=oId,subject=tId,
	id = postTTS_ir $id[tId],
	fiveSeq_id= fiveSeq_not_within_biprom$id[oId]
)

fiveSeq_within_postTTS = fiveSeq_not_within_biprom[ which( fiveSeq_not_within_biprom $id %in% minfo_withOL_postTTS$fiveSeq_id ), ]
fiveSeq_not_within_postTTS = fiveSeq_not_within_biprom[ which( !( fiveSeq_not_within_biprom $id %in% minfo_withOL_postTTS$fiveSeq_id) ), ]


## Antisense of gene

oltif_within_asgene = findOverlaps(fiveSeq_not_within_postTTS, gene_antisense_ir)
mfo = as.matrix(oltif_within_asgene)
oId=mfo[,1]
tId=mfo[,2]
minfo_withOL_asgene=data.frame(query=oId,subject=tId,
	gene_id = gene_antisense_ir$gene[tId],
	fiveSeq_id= fiveSeq_not_within_postTTS$id[oId]
)

fiveSeq_within_antisense = fiveSeq_not_within_postTTS[ which( fiveSeq_not_within_postTTS$id %in% minfo_withOL_asgene$fiveSeq_id ), ]
fiveSeq_not_within_antisense = fiveSeq_not_within_postTTS[ which( !( fiveSeq_not_within_postTTS$id %in% minfo_withOL_asgene$fiveSeq_id) ), ]



####### PLOT

### ONE PDF WITH POSITIONS

TSS_df = as.data.frame(fiveSeq_within_TSS)
UTR5_df = as.data.frame(fiveSeq_within_5UTR)
UTR3_df = as.data.frame(fiveSeq_within_3UTR)
exon_df = as.data.frame(fiveSeq_within_exon)
intron_df = as.data.frame(fiveSeq_within_intron)
biprom_df = as.data.frame(fiveSeq_within_biprom)
postTTS_df = as.data.frame(fiveSeq_within_postTTS)
as_df = as.data.frame(fiveSeq_within_antisense)
notAs_df = as.data.frame(fiveSeq_not_within_antisense)


pdf(outPdfPosition, width=15, height=15)
par(mfrow=c(ceiling(length(samples)/2),2))
#par(mfrow=c(2, ceiling(length(samples)/2)))

for(s in samples) {

	a=length(TSS_df[which(TSS_df[,grep(s, colnames(TSS_df))] >0),]$chr)
	b=length(UTR5_df[which(UTR5_df[,grep(s, colnames(UTR5_df))] >0),]$chr)
	c=length(UTR3_df[which(UTR3_df[,grep(s, colnames(UTR3_df))] >0),]$chr)
	d=length(exon_df[which(exon_df[,grep(s, colnames(exon_df))] >0),]$chr)
	e=length(intron_df[which(intron_df[,grep(s, colnames(intron_df))] >0),]$chr)
	f=length(biprom_df[which(biprom_df[,grep(s, colnames(biprom_df))] >0),]$chr)
	g=length(postTTS_df[which(postTTS_df[,grep(s, colnames(postTTS_df))] >0),]$chr)
	h=length(as_df[which(as_df[,grep(s, colnames(as_df))] >0),]$chr)
	i=length(notAs_df[which(notAs_df[,grep(s, colnames(notAs_df))] >0),]$chr)
	x = sum( c(a,b,c,d,e,f,g,h,i) )

	pie( 
		c(a,b,c,d,e,f,g,h,i),
		labels=c(
			paste( 'TSS \n(', round(100*a/x, digits=1), '%)', sep=''),
			paste( '5UTR \n(', round(100*b/x, digits=1), '%)', sep=''),
			paste( '3UTR \n(', round(100*c/x, digits=1), '%)', sep=''),
			paste( 'exon \n(', round(100*d/x, digits=1), '%)', sep=''),
			paste( 'intron \n(', round(100*e/x, digits=1), '%)', sep=''),
			paste( 'bi-prom (', round(100*f/x, digits=1), '%)', sep=''),
			paste( 'post-TTS (', round(100*g/x, digits=1), '%)', sep=''),
			paste( 'antisense (', round(100*h/x, digits=1), '%)', sep=''),
			paste( 'intergenic \n(', round(100*i/x, digits=1), '%)', sep='')
	),
	col=brewer.pal(9, "Set3"), # [1:length(samples)],
	main=paste("Assignment of 5Seq positions\n(", s, ")", sep='')
	)
}

dev.off()


### ANOTHER PDF WITH COUNTS


pdf(outPdfCount, width=15, height=15)
par(mfrow=c(ceiling(length(samples)/2),2))
#par(mfrow=c(2, ceiling(length(samples)/2)))

for(s in samples) {

	a=sum(TSS_df[which(TSS_df[,grep(s, colnames(TSS_df))] >0),grep(s, colnames(TSS_df))])
	b=sum(UTR5_df[which(UTR5_df[,grep(s, colnames(UTR5_df))] >0),grep(s, colnames(UTR5_df))])
	c=sum(UTR3_df[which(UTR3_df[,grep(s, colnames(UTR3_df))] >0),grep(s, colnames(UTR3_df))])
	d=sum(exon_df[which(exon_df[,grep(s, colnames(exon_df))] >0),grep(s, colnames(exon_df))])
	e=sum(intron_df[which(intron_df[,grep(s, colnames(intron_df))] >0),grep(s, colnames(intron_df))])
	f=sum(biprom_df[which(biprom_df[,grep(s, colnames(biprom_df))] >0),grep(s, colnames(biprom_df))])
	g=sum(postTTS_df[which(postTTS_df[,grep(s, colnames(postTTS_df))] >0),grep(s, colnames(postTTS_df))])
	h=sum(as_df[which(as_df[,grep(s, colnames(as_df))] >0),grep(s, colnames(as_df))])
	i=sum(notAs_df[which(notAs_df[,grep(s, colnames(notAs_df))] >0),grep(s, colnames(notAs_df))])
	x = sum( c(a,b,c,d,e,f,g,h,i) )

	pie( 
		c(a,b,c,d,e,f,g,h,i),
		labels=c(
			paste( 'TSS \n(', round(100*a/x, digits=1), '%)', sep=''),
			paste( '5UTR \n(', round(100*b/x, digits=1), '%)', sep=''),
			paste( '3UTR \n(', round(100*c/x, digits=1), '%)', sep=''),
			paste( 'exon \n(', round(100*d/x, digits=1), '%)', sep=''),
			paste( 'intron \n(', round(100*e/x, digits=1), '%)', sep=''),
			paste( 'bi-prom (', round(100*f/x, digits=1), '%)', sep=''),
			paste( 'post-TTS (', round(100*g/x, digits=1), '%)', sep=''),
			paste( 'antisense (', round(100*h/x, digits=1), '%)', sep=''),
			paste( 'intergenic \n(', round(100*i/x, digits=1), '%)', sep='')
	),
	col=brewer.pal(9, "Set3"), # [1:length(samples)],
	main=paste("Assignment of 5Seq counts\n(", s, ")", sep='')
	)
}

dev.off()


