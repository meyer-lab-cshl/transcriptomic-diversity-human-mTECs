###############################
## libraries and functions ####
###############################

library(AnnotationHub)
library(data.table)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9)
library(dplyr)

preprocess_feature_table <- function(datx, daty) {
    # merge dataframes to find non-overlapping positions
    dat <- merge(datx, daty, by=c("geneID","Position"), all.x = TRUE )
    dat <- subset(dat, select = -c(ReadCount.x, PosFromAnno.x, Class.x,
                                   ReadCount.y, PosFromAnno.y, Class.y))

    # define response value: In other data source? yes/no
    dat$response<- "blub"
    dat[is.na(dat$BarcodeCount.y),]$response<- 'no'
    dat[!is.na(dat$BarcodeCount.y),]$response <- 'yes'

    # split into genes and intergenic regions
    dat_inter <- dat[grepl("intergenic", dat$geneID),]
    dat <- dat[!grepl("intergenic", dat$geneID),]

    # for genes annotate chromosome and strand with anno_genes object
    dat$chrom <- with(dat,
                      as.character(seqnames(anno_genes[as.character(dat$geneID)])))
    dat$strand <- with(dat, as.character(strand(anno_genes[dat$geneID])))

    # for intergenic regions annotate chromosome from name and strand as plus
    dat_inter$chrom <- gsub(".*_(chr.{1,2})_.*", "\\1", dat_inter$geneID)
    dat_inter$strand <- gsub(".*_.*_(.)_.*", "\\1", dat_inter$geneID)

    # merge dataframes again
    dat <- rbind(dat, dat_inter)

    # define end
    dat$end <- as.integer(dat$Position +1)
    dat <- dat[c(1,6,7,2,8,5,3)]
    return(dat)
}

generate_features <- function(anno_exons, anno_genes, anno_transcripts) {
    features <- list()
    # Define GRanges
    ## TSS +-10
    features$TSS10_genes <- promoters(anno_genes, upstream=10, downstream=10)

    ## TSS anno_transcriptss +-10
    features$TSS10_anno_transcriptss <- promoters(anno_transcripts, upstream=10,
                                            downstream=10)

    ## TSS +-100
    features$TSS100_genes <- promoters(anno_genes, upstream=100, downstream=100)

    ## TSS transcrips +-10
    features$TSS100_transcripts <- promoters(anno_transcripts, upstream=100,
                                             downstream=100)

    ## First Exon
    # exon is genomicranges list object, for each element of list take either
    # first or last entry, depending on strand
    first_exon <- lapply(anno_exons, function(x) {
                         if( unique(strand(x)) == '+') x[1] else x[length(x)]
                  })
    first_exon <- do.call(GRangesList, first_exon)
    features$first_exon <- unlist(first_exon)

    ## Other Exons
    # same as first exon, except all exon ranges are taken except for first or
    # last entry, depending on strand
    other_exon <- lapply(anno_exon, function(x) {
                        if( unique(strand(x)) == '+') x[-1] else x[-length(x)]
                      })
    other_exon <- do.call(GRangesList, other_exon)
    features$other_exon <- unlist(other_exon)

    ## Introns
    # later a different more correct definition was used
    all_exons <- unlist(anno_exon)
    features$introns <- GenomicRanges::setdiff(anno_genes, all_exons)

    ## everything from here on has to be strand specific
    ## TTS +-100
    plus <- anno_genes[strand(anno_genes) == '+']
    minus <- anno_genes[strand(anno_genes) == '-']
    plus_tts <- GRanges(seqnames=seqnames(plus),
                        ranges=IRanges(start=(end(plus)-100),
                                       end=(end(plus)+100)),
                        strand=strand(plus))
    minus_tts <- GRanges(seqnames=seqnames(minus),
                         ranges=IRanges(start=(start(minus)-100),
                                        end=(start(minus)+100)),
                         strand = strand(minus))
    features$TTS100_genes <- c(plus_tts, minus_tts)

    ## TTS +200
    start(plus_tts) <- start(plus_tts) + 100
    end(plus_tts) <- end(plus_tts) + 100
    start(minus_tts) <- start(minus_tts) - 100
    end(minus_tts) <- end(minus_tts) - 100
    features$TTS200_genes <- c(plus_tts, minus_tts)

    ## TTS transcripts +-100
    plus_trans <- anno_transcripts[strand(anno_transcripts) == '+']
    minus_trans <- anno_transcripts[strand(anno_transcripts) == '-']
    plus_trans_tts <- GRanges(seqnames=seqnames(plus_trans),
                              ranges=IRanges(start=(end(plus_trans)-100),
                                             end=(end(plus_trans)+100)),
                              strand = strand(plus_trans))
    minus_trans_tts <- GRanges(seqnames=seqnames(minus_trans),
                               ranges=IRanges(start=(start(minus_trans)-100),
                                              end=(start(minus_trans)+100)),
                               strand=strand(minus_trans))
    features$TTS100_trans <- c(plus_trans_tts, minus_trans_tts)

    ## TTS transcripts +200
    start(plus_trans_tts) <- start(plus_trans_tts) + 100
    end(plus_trans_tts) <- end(plus_trans_tts) + 100
    start(minus_trans_tts) <- start(minus_trans_tts) - 100
    end(minus_trans_tts) <- end(minus_trans_tts) - 100
    features$TTS200_trans <- c(plus_trans_tts, minus_trans_tts)

    ## Downstream
    plus <- anno_genes[strand(anno_genes) == '+']
    minus <- anno_genes[strand(anno_genes) == '-']
    start(plus) <- end(plus)
    end(plus) <- end(plus) + 1000
    end(minus) <- start(minus)
    start(minus) <- start(minus) - 1000
    features$downstream_genes <- c(plus, minus)

    ## Upstream
    plus <- anno_genes[strand(anno_genes) == '+']
    minus <- anno_genes[strand(anno_genes) == '-']
    end(plus) <- start(plus)
    start(plus) <- start(plus) - 1000
    start(minus) <- end(minus)
    end(minus) <- end(minus) + 1000
    features$upstream_genes <- c(plus, minus)

    ## Antisense
    # later this object was extended by 1000
    plus <- anno_genes[strand(anno_genes) == '+']
    minus <- anno_genes[strand(anno_genes) == '-']
    strand(plus) <- '-'
    strand(minus) <- '+'
    features$antisense <- c(plus, minus)

    return(features)
}

collect_features <- function(total, feature_list) {
    ## Rank genes by expression level ####
    genes <- total[,c("geneID", "start", "BarcodeCount.x")]

    # sum up count from all positions of each gene
    genes <- genes[,lapply(.SD, sum), by=.(geneID), .SDcols=c("BarcodeCount.x")]
    setDF(genes)

    # define a rank for each gene based on overall expression level
    genes$rank <- rank(genes$BarcodeCount.x, ties.method = "average")
    total$expr_rank_gene <- with(total, genes$rank[geneID])

    ## Expressionlevel relativ to gene and to total for each peak ####
    total <- total %>%
        group_by(geneID) %>%
        mutate(sumCount=sum(BarcodeCount.x)) %>%
        mutate(relCount=BarcodeCount.x/sumCount) %>%
        group_by() %>%
        mutate(relAllCount=as.numeric(BarcodeCount.x/sum(BarcodeCount.x))) %>%
        dplyr::select(-sumCount)

    ## Annotate features ####
    ranges <- makeGRangesFromDataFrame(total, keep.extra.columns = TRUE)

    total$TSS_10_genes <- overlapsAny(ranges, TSS10_genes)
    total$TSS_10_trans <- overlapsAny(ranges, TSS10_anno_transcriptss)
    total$TSS_100_genes <- overlapsAny(ranges, TSS100_genes)
    total$TSS_100_trans <- overlapsAny(ranges, TSS100_anno_transcriptss)
    total$first_exon <- overlapsAny(ranges, first_exon)
    total$other_exon <- overlapsAny(ranges, other_exon)
    total$intron <- overlapsAny(ranges, introns)
    total$TTS_100_genes <- overlapsAny(ranges, TTS100_genes)
    total$TTS_200_genes <- overlapsAny(ranges, TTS200_genes)
    total$TTS_100_trans <- overlapsAny(ranges, TTS100_trans)
    total$TTS_200_trans <- overlapsAny(ranges, TTS200_trans)
    total$downstream_gene <- overlapsAny(ranges, downstream_genes)
    total$upstream_gene <- overlapsAny(ranges, upstream_genes)
    total$antisense <- overlapsAny(ranges, antisense)

    ## Base content ####
    # getSeq to get actual sequence from genome object containing genomic
    # sequence, reverse for minus strand
    total_plus <- total[total$strand == '+',]
    total_minus <- total[total$strand == '-',]
    total_plus$Basesafter <- getSeq(genome, total_plus$chrom,
                                    (total_plus$start + 1),
                                    (total_plus$start + 2),
                                    as.character=TRUE)
    total_minus$Basesafter <- reverse(getSeq(genome, total_minus$chrom,
                                            (total_minus$start - 2),
                                            (total_minus$start - 1),
                                            as.character=TRUE))
    total <- rbind(total_plus,total_minus)

    ## Base content 10bp window
    total$Freq10A <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-5),
                                            (total$start+5)), "A")
    total$Freq10C <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-5),
                                            (total$start+5)), "C")
    total$Freq10G <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-5),
                                            (total$start+5)), "G")
    total$Freq10T <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-5),
                                            (total$start+5)), "T")

    ## Base content 50bp window
    total$Freq50A <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-25),
                                            (total$start+25)), "A")
    total$Freq50C <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-25),
                                            (total$start+25)), "C")
    total$Freq50G <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-25),
                                            (total$start+25)), "G")
    total$Freq50T <- letterFrequency(getSeq(genome, total$chrom,
                                            (total$start-25),
                                            (total_int$start+25)), "T")
    return(total)
}

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-m", "--m5pseqdir"), action="store", dest="mousedir",
               type="character", help="Path to 5Pseq mouse isoforms directory
               [default: %default].", default=NULL),
    make_option(c("-f", "--fantomdir"), action="store", dest="fantomdir",
               type="character", help="Path to Fantom5 mouse isoforms directory
               [default: %default].", default=NULL),
    make_option(c("-o", "--odir"), action="store", dest="odir",
               type="character", help="Path to output directory
               [default: %default].", default=NULL),
    make_option(c("-s", "--suffix"), action="store", dest="suffix",
               type="character", help="common suffix of files to read
               [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$fantom <- "~/data/tss/combined/tss/isoforms/all_mouse_ES_46C_fantoms.isoforms.csv"
    args$mouse <- "~/data/tss/combined/tss/isoforms/all_mouse_mESC.isoforms.csv"
    args$odir <- "~/data/tss/combined/tss/isoforms"
    args$suffix <- ".isoforms.csv"
}

# define mouse chromosome ids
chr <- paste("chr", c(1:22, "M", "X", "Y"), sep="")

## Load annotations ####
hub <- AnnotationHub()
# |Organism: Mus musculus
# |taxonomy_id: 10090
# |genome_build: GRCm38
# |DBSCHEMAVERSION: 2.1
# | No. of genes: 56289.
# | No. of anno_transcriptss: 144726.
anno <- hub[["AH78811"]]
seqlevelsStyle(anno) <- "UCSC"
chr <- paste ("chr", c(1:19, "M", "X", "Y"), sep="")

# get gene ids, expand range by 100
anno_genes <- genes(anno)
anno_transcripts <- anno_transcripts(anno)
anno_exon <- exonsBy(anno, by='gene')
anno_genes <- subset(anno_genes, seqnames %in% chr)
anno_transcripts <- subset(anno_transcripts, seqnames %in% chr)

genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9

# read TSS data
m5pseq <- fread(args$mouse, sep=",", header = TRUE, stringsAsFactors=TRUE,
                  data.table=FALSE)
fantom <- fread(args$fantom, sep=",", header = TRUE, stringsAsFactors=TRUE,
                  data.table=FALSE)

################
## analysis ####
################

## Generate feature table for annotation ####
features <- generate_features(anno_exons, anno_genes, anno_transcripts)

## Pre-process tss data ####
total_m5pseq <- preprocess_feature_table(m5pseq, fantom)
total_m5pseq <- dplyr::rename(total_m5pseq, start=Position, Fantom5=response)

total_fantom <- preprocess_feature_table(fantom, internal)
total_fantom <- dplyr::rename(total_fantom, start=Position, Internal=response)

## Annotate tss data with features ####
features_fantom <- collect_features(total_fantom, features)
features_m5pseq <- collect_features(total_m5pseq, features)

write.table(features_m5pseq,
            file=file.path(args$odir, "features_5pseq_mm9.csv"),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(features_fantom,
            file=file.path(args$odir, "features_fantom5_mm9.csv"),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
