# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      Preparing various annotation objects for assignment.
#
# cat /g/steinmetz/project/mTEC_Seq/src/prepareAssignmentObjects.R | R --slave --args
#	gtf='/g/steinmetz/genome/Homo_sapiens/37.68/annotation/gtf/Homo_sapiens.GRCh37.68.chrOnly.gtf'
#	outfile='/g/steinmetz/project/mTEC_Seq/annotation/assignment_IRanges.rda'
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


####### FIXED PARAMETERS

genomic=c( '1','10','11','12','13','14','15','16','17','18','19','2','3','4','5','6','7','8','9','20','21','22','X','Y' )



####### PROCESS THE GTF FILE TO YIELD VARIOUS ANNOTATION OBJECTS

gtf = read.delim(gtf,header=F)

exons = gtf[which(gtf$V3=='exon'), ]
exons = exons[which(as.character(exons$V1) %in% genomic),]


## ## ## ## All unique TSS:

# For each transcript, get start position
txids = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'transcript_id ', strsplit( x, '; ' )[[1]] )], 'transcript_id ' )[[1]][2] ) )


a=cbind(exons, txids)

plus=a[which(a[,7]=='+'),]
minus=a[which(a[,7]=='-'),]
plus_starts = by(plus[,4], as.character(plus$txids), min)
minus_starts = by(minus[,5], as.character(minus$txids), max)
plus_chrs = by(plus[,1], as.character(plus$txids), function(x) { as.character(x[1]) } )
minus_chrs = by(minus[,1], as.character(minus$txids), function(x) { as.character(x[1]) } )
plus_strands = by(plus[,7], as.character(plus$txids), function(x) { as.character(x[1]) } )
minus_strands = by(minus[,7], as.character(minus$txids), function(x) { as.character(x[1]) } )

t = as.data.frame(matrix(nrow=length(plus_starts) + length(minus_starts), ncol=0))
t$chr = c(plus_chrs, minus_chrs)
t$strand = c(plus_strands, minus_strands)
t$TSS = c(plus_starts, minus_starts)
t$txid = c(names(plus_strands), names(minus_strands))
t2=unique(t[,1:3])


# IR:
TSS_ir = RangedData(
	IRanges(
	start= t2$TSS,
	end= t2$TSS),
	space=paste( t2$chr,'.', t2$strand,sep=''),
	chr= t2$chr,
	strand= t2$strand,
	id=paste(t2$chr, t2$strand, t2$TSS ,sep='_')
)



## ## ## ## Exons

x = unique( exons[, c(1,4,5,7)] )

exon_ir = RangedData(
	IRanges(
	start= x[,2],
	end= x[,3]),
	space=paste( x[,1],'.', x[,4],sep=''),
	chr= x[,1],
	strand= x[,4]
	#id=paste(x[,1], [,], t2$TSS ,sep='_')
)


## ## ## ## Full gene (intronic will be those assigning to gene but NOT to exon)

geneids = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'gene_id ', strsplit( x, '; ' )[[1]] )], 'gene_id ' )[[1]][2] ) )

a=cbind(exons, geneids)

starts = by(a[,4], as.character(a$geneids), min )
ends = by(a[,5], as.character(a$geneids), max )
chrs = by(a[,1], as.character(a$geneids), function(x) { as.character(x[1]) } )
strands = by(a[,7], as.character(a$geneids), function(x) { as.character(x[1]) } )

t = as.data.frame(matrix(nrow=length(starts), ncol=0))
t$chr = as.character(chrs)
t$strand = as.character(strands)
t$start = as.integer(starts)
t$end = as.integer(ends)
t$geneid = names(ends)

gene_ir = RangedData(
	IRanges(
	start= t$start,
	end= t$end),
	space=paste( t$chr,'.', t$strand,sep=''),
	chr= t$chr,
	strand= t$strand,
	gene=t$geneid
)


## ## ## ## Antisense (of full gene)

gene_antisense_ir = RangedData(
	IRanges(
	start= t$start,
	end= t$end),
	space=paste( t$chr,'.', ifelse( t$strand=='+', '-','+' ),sep=''),
	chr= t$chr,
	strand= t$strand,
	gene=t$geneid
)



## ## ## ## 5UTR and 3UTR annotation (when taking not within +/-x TSS, this can be entire 5UTR)

# For each transcript, get all exons and coding region
# Select exons "before" coding region. 
# From the exon with start codon, keep positions before coding region.
# Make a table and add to IR (row per exons coordinates): chr strand start end txid (hasStartCodon)

txids = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'transcript_id ', strsplit( x, '; ' )[[1]] )], 'transcript_id ' )[[1]][2] ) )
geneids = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'gene_id ', strsplit( x, '; ' )[[1]] )], 'gene_id ' )[[1]][2] ) )
exonN = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'exon_number ', strsplit( x, '; ' )[[1]] )], 'exon_number ' )[[1]][2] ) )
a=cbind(exons[,c(1,4,5,7)], txids, geneids, exonN)

## CDS table: txid startCodonPos stopCodonPos startCodonExon stopCodonExon
start_codons = gtf[which(gtf$V3=='start_codon'),]
stop_codons = gtf[which(gtf$V3=='stop_codon'),]

start_codons_txids = unlist( lapply(as.character( start_codons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'transcript_id ', strsplit( x, '; ' )[[1]] )], 'transcript_id ' )[[1]][2] ) )
stop_codons_txids = unlist( lapply(as.character( stop_codons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'transcript_id ', strsplit( x, '; ' )[[1]] )], 'transcript_id ' )[[1]][2] ) )

start_codons_exon = unlist( lapply(as.character( start_codons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'exon_number ', strsplit( x, '; ' )[[1]] )], 'exon_number ' )[[1]][2] ) )
stop_codons_exon = unlist( lapply(as.character( stop_codons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'exon_number ', strsplit( x, '; ' )[[1]] )], 'exon_number ' )[[1]][2] ) )

start_codons_geneids = unlist( lapply(as.character( start_codons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'gene_id ', strsplit( x, '; ' )[[1]] )], 'gene_id ' )[[1]][2] ) )
stop_codons_geneids = unlist( lapply(as.character( stop_codons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'gene_id ', strsplit( x, '; ' )[[1]] )], 'gene_id ' )[[1]][2] ) )

startCodons = cbind(start_codons[,c(1,7,4,5)], start_codons_geneids, start_codons_txids, start_codons_exon)
stopCodons = cbind(stop_codons[,c(1,7,4,5)], stop_codons_geneids, stop_codons_txids, stop_codons_exon)

colnames(startCodons) = c('chr', 'strand', 'start', 'end', 'geneid', 'txid', 'exonWithStartCodon')
colnames(stopCodons) = c('chr', 'strand', 'start', 'end', 'geneid', 'txid', 'exonWithStopCodon')

## Merge all the tables:
cds_coords = merge(startCodons, stopCodons, by.x='txid', by.y='txid')

# Just to check:
#table(cds_coords$chr.x == cds_coords$chr.y)
#table(cds_coords$strand.x == cds_coords$strand.y)

b = merge(a, cds_coords, by.x='txids', by.y='txid')

txt = as.data.frame(matrix(nrow=length(b[,1]), ncol=0))
txt$chr = as.character( b$V1 )
txt$strand = as.character( b$V7 )
txt$start = b$V4
txt$end = b$V5
txt$exon = as.character( b$exonN )
txt$startCodonStart = b$start.x
txt$startCodonEnd = b$end.x
txt$stopCodonStart = b$start.y
txt$stopCodonEnd = b$end.y
txt$txid = as.character( b$txids )
txt$exonWithStartCodon = as.integer( as.character( b$exonWithStartCodon ) )
txt$exonWithStopCodon = as.integer( as.character( b$exonWithStopCodon ) )

#table( as.integer(as.character(txt$exonWithStartCodon)) <= as.integer(as.character(txt$exonWithStopCodon)) )
#  TRUE 
#684055

# For each transcript, select exons (by exon number) that are smaller or equal to exonWithStartCodon (--> 5UTR exons).
# 'Chop up' the start codon containing exon, and make the 5UTR table.

txtsplit = split(txt, as.character(txt$txid))

# For TX with multiple start codons, use the first one.
nStartCodons = lapply(txtsplit, function(x) {
	return( length(unique(x$exonWithStartCodon)) )
})

zz = do.call(rbind, nStartCodons)
#table(zz[,1])
#    1     2 
#64630   170


utr5info = lapply(txtsplit, function(x) {
	# Get 5UTR exons	
	x = unique(x[,c(1:7, 10:11)])
	UTR5exons = x[ which( as.integer(as.character(x[,'exon'])) <= as.integer(as.character(x[,'exonWithStartCodon']))), ]
	# For the exon WITH start codon
	codonexon = UTR5exons[which(UTR5exons$exon==UTR5exons$exonWithStartCodon),]
	# If multiple start codons in this transcript, choose based on the earliest one
	if( length(unique(codonexon$exonWithStartCodon))>1 ) {
		#
		if( codonexon$strand[1] == '+' ) {
			codonexon = codonexon[ which( codonexon$startCodonStart == min(codonexon$startCodonStart) ), ][1,]
		} else {
			codonexon = codonexon[ which( codonexon$startCodonEnd == max(codonexon$startCodonEnd) ), ][1,]
		}
	}

	# Get positions UNTIL start codon (mind the strand)
	# Always choosing the first instance because
	if(codonexon$strand == '+') {
		codonexon$end = codonexon$startCodonStart - 1
	} else {
		codonexon$start = codonexon$startCodonEnd + 1
	}

	# The rest of the exons
	noncodonexon = UTR5exons[which(UTR5exons$exon!=UTR5exons$exonWithStartCodon),]
	# Multi-start-codon-transcripts
	if( length(unique(codonexon$exonWithStartCodon))>1 ) {
		#
		if( codonexon$strand[1] == '+' ) {
			additionalexon = codonexon[ which( codonexon$startCodonStart != min(codonexon$startCodonStart) ), ]
		} else {
			additionalexon = codonexon[ which( codonexon$startCodonEnd != max(codonexon$startCodonEnd) ), ]
		}
		noncodonexon = rbind(noncodonexon, additionalexon)
	}

	# rbind and return
	return(rbind(noncodonexon, codonexon))
})


utr5infoT = do.call(rbind,utr5info)
utr5infoT2 = utr5infoT[order(utr5infoT$chr, utr5infoT$strand, utr5infoT$start, utr5infoT$end, utr5infoT$txid), ]


# For TX with multiple start codons, use the first one.
nStopCodons = lapply(txtsplit, function(x) {
	return( length(unique(x$exonWithStopCodon)) )
})

yy = do.call(rbind, nStopCodons)
# table(yy[,1])

utr3info = lapply(txtsplit, function(x) {
	x = unique(x[,c(1:5, 8:10, 12)])

	# Get 3UTR exons	
	UTR3exons = x[ which( as.integer(as.character(x[,'exon'])) >= as.integer(as.character(x[,'exonWithStopCodon']))), ]
	# For the exon WITH stop codon
	codonexon = UTR3exons[which(UTR3exons$exon==UTR3exons$exonWithStopCodon),]

	# If multiple stop codons in this transcript, choose based on the last one
	if( length(unique(codonexon$exonWithStopCodon))>1 ) {
		#
		if( codonexon$strand[1] == '+' ) {
			codonexon = codonexon[ which( codonexon$stopCodonEnd == max(codonexon$stopCodonEnd) ), ][1,]
		} else {
			codonexon = codonexon[ which( codonexon$stopCodonStart == min(codonexon$stopCodonStart) ), ][1,]
		}
	}

	# Get positions UNTIL stop codon (mind the strand)
	# 
	if(codonexon$strand == '+') {
		codonexon$start = codonexon$stopCodonEnd + 1
	} else {
		codonexon$end = codonexon$stopCodonStart - 1
	}

	# The rest of the exons
	noncodonexon = UTR3exons[which(UTR3exons$exon!=UTR3exons$exonWithStopCodon),]	

	# Multi-stop-codon-transcripts
	if( length(unique(codonexon$exonWithStopCodon))>1 ) {
		#
		if( codonexon$strand[1] == '+' ) {
			additionalexon = codonexon[ which( codonexon$stopCodonEnd == max(codonexon$stopCodonEnd) ), ][1,]
		} else {
			additionalexon = codonexon[ which( codonexon$stopCodonStart == min(codonexon$stopCodonStart) ), ][1,]
		}
		noncodonexon = rbind(noncodonexon, additionalexon)
	}


	# rbind and return
	return(rbind(noncodonexon, codonexon))
})

utr3infoT = do.call(rbind,utr3info)
utr3infoT2 = utr3infoT[order(utr3infoT$chr, utr3infoT$strand, utr3infoT$start, utr3infoT$end, utr3infoT$txid), ]


## Build an IR objects

fiveUTR_ir = RangedData(
	IRanges(
	start= utr5infoT2$start,
	end= utr5infoT2$end),
	space=paste( utr5infoT2$chr,'.', utr5infoT2$strand,sep=''),
	chr= utr5infoT2$chr,
	strand= utr5infoT2$strand,
	txid=utr5infoT2$txid,
	exon=utr5infoT2$exon,
	exonWithStartCodon =utr5infoT2$exonWithStartCodon,
	startCodonStart=utr5infoT2$startCodonStart,
	startCodonEnd=utr5infoT2$startCodonEnd
)

threeUTR_ir = RangedData(
	IRanges(
	start= utr3infoT2$start,
	end= utr3infoT2$end),
	space=paste( utr3infoT2$chr,'.', utr3infoT2$strand,sep=''),
	chr= utr3infoT2$chr,
	strand= utr3infoT2$strand,
	txid=utr3infoT2$txid,
	exon=utr3infoT2$exon,
	exonWithStopCodon =utr3infoT2$exonWithStopCodon,
	stopCodonStart=utr3infoT2$stopCodonStart,
	stopCodonEnd=utr3infoT2$stopCodonEnd
)




## ## ## ## Bidirectional promoter (antisense from TSS, within 1kb)


# For each transcript, get start position
txids = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'transcript_id ', strsplit( x, '; ' )[[1]] )], 'transcript_id ' )[[1]][2] ) )

a=cbind(exons, txids)

plus=a[which(a[,7]=='+'),]
minus=a[which(a[,7]=='-'),]
plus_starts = by(plus[,4], as.character(plus$txids), min)
minus_starts = by(minus[,5], as.character(minus$txids), max)
plus_chrs = by(plus[,1], as.character(plus$txids), function(x) { as.character(x[1]) } )
minus_chrs = by(minus[,1], as.character(minus$txids), function(x) { as.character(x[1]) } )
plus_strands = by(plus[,7], as.character(plus$txids), function(x) { as.character(x[1]) } )
minus_strands = by(minus[,7], as.character(minus$txids), function(x) { as.character(x[1]) } )

t = as.data.frame(matrix(nrow=length(plus_starts) + length(minus_starts), ncol=0))
t$chr = c(plus_chrs, minus_chrs)
t$strand = c(plus_strands, minus_strands)
t$TSS = c(plus_starts, minus_starts)
t$txid = c(names(plus_strands), names(minus_strands))

t2=unique(t[,1:3])

# IR:
extension=1000
biprom_ir = RangedData(
	IRanges(
	start= ifelse( t2$strand=='+', t2$TSS - extension, t2$TSS + 1),
	end= ifelse( t2$strand=='+', t2$TSS - 1, t2$TSS + extension )),
	space=paste( t2$chr,'.', ifelse( t2$strand=='+', '-', '+' ),sep=''),
	chr= t2$chr,
	strand= t2$strand,
	id=paste(t2$chr, t2$strand, t2$TSS ,sep='_')
)


## ## ## ## Soon after gene (sense strand, within 1kb)

# For each transcript, get END position
txids = unlist( lapply(as.character( exons[,9] ), function(x) strsplit( strsplit( x, '; ' )[[1]][ grep( 'transcript_id ', strsplit( x, '; ' )[[1]] )], 'transcript_id ' )[[1]][2] ) )
a=cbind(exons, txids)

plus=a[which(a[,7]=='+'),]
minus=a[which(a[,7]=='-'),]
plus_ends = by(plus[,5], as.character(plus$txids), max)
minus_ends = by(minus[,4], as.character(minus$txids), min)
plus_chrs = by(plus[,1], as.character(plus$txids), function(x) { as.character(x[1]) } )
minus_chrs = by(minus[,1], as.character(minus$txids), function(x) { as.character(x[1]) } )
plus_strands = by(plus[,7], as.character(plus$txids), function(x) { as.character(x[1]) } )
minus_strands = by(minus[,7], as.character(minus$txids), function(x) { as.character(x[1]) } )

t = as.data.frame(matrix(nrow=length(plus_starts) + length(minus_starts), ncol=0))
t$chr = c(plus_chrs, minus_chrs)
t$strand = c(plus_strands, minus_strands)
t$TTS = c(plus_ends, minus_ends)
t$txid = c(names(plus_strands), names(minus_strands))

t2=unique(t[,1:3])


# Make a table, unique.

# IR:
extension=1000
postTTS_ir = RangedData(
	IRanges(
	start= ifelse( t2$strand=='+', t2$TTS + 1, t2$TTS - extension ),
	end= ifelse( t2$strand=='+', t2$TTS + extension, t2$TTS - 1 ) ),
	space=paste( t2$chr,'.', t2$strand,sep=''),
	chr= t2$chr,
	strand= t2$strand,
	id=paste(t2$chr, t2$strand, t2$TTS ,sep='_')
)


save(TSS_ir, exon_ir, gene_ir, gene_antisense_ir, fiveUTR_ir, threeUTR_ir, biprom_ir, postTTS_ir, file=outfile)


