###############################
## libraries and functions ####
###############################

library(AnnotationHub)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)

preprocess_feature_table <- function(datx, daty, anno_genes) {
    # merge dataframes to find non-overlapping positions
    # define response value: In other data source? yes/no
    dat <- left_join(datx, daty, by=c("geneID", "Position")) %>%
        arrange(geneID, Position) %>%
        dplyr::select(geneID, Position, BarcodeCount.x, BarcodeCount.y) %>%
        mutate(response=case_when(is.na(BarcodeCount.y) ~ "no",
                                        TRUE ~ "yes"))

    # split into genes and intergenic regions
    dat_inter <- dat[grepl("intergenic", dat$geneID),]
    dat <- dat[!grepl("intergenic", dat$geneID),]

    # for genes annotate chromosome and strand with anno_genes object
    dat <- dat %>%
        mutate(chrom=as.character(seqnames(anno_genes[geneID])),
               strand=as.character(strand(anno_genes[geneID])))

    # for intergenic regions annotate chromosome from name and strand as plus
    dat_inter <- dat_inter %>%
        mutate(chrom=gsub(".*_(chr.{1,2})_.*", "\\1", geneID),
               strand=gsub(".*_.*_(.)_.*", "\\1", geneID))

    # merge dataframes again
    dat <- bind_rows(dat, dat_inter)

    # define end
    # re-order columns
    dat <- dat %>%
        mutate(end=Position + 1) %>%
        dplyr::select(geneID, chrom, strand, Position, end, response,
                      BarcodeCount.x)
    return(dat)
}

generate_features <- function(anno_exons, anno_genes, anno_transcripts) {
    features <- list()
    # Define GRanges
    ## TSS +-10
    features$TSS10_genes <- promoters(anno_genes, upstream=10, downstream=10)

    ## TSS anno_transcriptss +-10
    features$TSS10_transcripts <- promoters(anno_transcripts, upstream=10,
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

collect_features <- function(total, features, genome) {
    ## Rank genes by expression level ####
    # sum up count from all positions of each gene
    # and define a rank for each gene based on overall expression level
    genes <- total %>%
        dplyr::select(geneID, start, BarcodeCount) %>%
        group_by(geneID) %>%
        summarise(counts=sum(BarcodeCount)) %>%
        ungroup %>%
        mutate(expr_rank_gene = rank(counts, ties.method = "average"))

    ## Expression level relativ to gene and to total for each peak ####
    total <- total %>%
        group_by(geneID) %>%
        mutate(sumCount=sum(BarcodeCount)) %>%
        mutate(relCount=BarcodeCount/sumCount) %>%
        ungroup %>%
        mutate(relAllCount=as.numeric(BarcodeCount/sum(BarcodeCount))) %>%
        left_join(genes, by="geneID") %>%
        dplyr::select(-sumCount, counts)

    ## Annotate features ####
    ranges <- makeGRangesFromDataFrame(total, keep.extra.columns = TRUE)

    total$TSS_10_genes <- overlapsAny(ranges, features$TSS10_genes)
    total$TSS_10_trans <- overlapsAny(ranges, features$TSS10_transcripts)
    total$TSS_100_genes <- overlapsAny(ranges, features$TSS100_genes)
    total$TSS_100_trans <- overlapsAny(ranges, features$TSS100_transcripts)
    total$first_exon <- overlapsAny(ranges, features$first_exon)
    total$other_exon <- overlapsAny(ranges, features$other_exon)
    total$intron <- overlapsAny(ranges, features$introns)
    total$TTS_100_genes <- overlapsAny(ranges, features$TTS100_genes)
    total$TTS_200_genes <- overlapsAny(ranges, features$TTS200_genes)
    total$TTS_100_trans <- overlapsAny(ranges, features$TTS100_trans)
    total$TTS_200_trans <- overlapsAny(ranges, features$TTS200_trans)
    total$downstream_gene <- overlapsAny(ranges, features$downstream_genes)
    total$upstream_gene <- overlapsAny(ranges, features$upstream_genes)
    total$antisense <- overlapsAny(ranges, features$antisense)

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
    total <- rbind(total_plus, total_minus)

    ## Base content 10bp window
    start_m5 <- total$start - 5
    start_m5[start_m5 < 1] <- 1
    total$Freq10A <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m5,
                                            (total$start+5)), "A")
    total$Freq10C <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m5,
                                            (total$start+5)), "C")
    total$Freq10G <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m5,
                                            (total$start+5)), "G")
    total$Freq10T <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m5,
                                            (total$start+5)), "T")

    ## Base content 50bp window
    start_m25 <- total$start - 25
    start_m25[start_m25 < 1] <- 1
    total$Freq50A <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m25,
                                            (total$start+25)), "A")
    total$Freq50C <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m25,
                                            (total$start+25)), "C")
    total$Freq50G <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m25,
                                            (total$start+25)), "G")
    total$Freq50T <- letterFrequency(getSeq(genome, total$chrom,
                                            start_m25,
                                            (total$start+25)), "T")
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

extendRange <- function(annotation, offset=100) {
    tmp <- tibble(chrom=as.character(seqnames(annotation)),
                  gene=annotation$gene_id,
           end_range=end(annotation),
           start_range=start(annotation))
    seqlen <- seqlengths(annotation) %>%
        enframe(name="chrom", value="length")
    tmp <- tmp %>%
        left_join(seqlen, by="chrom") %>%
        dplyr::filter(end_range + offset > length | start_range - offset < 0)

    annotation[!annotation$gene_id %in% tmp$gene,] <-
        annotation[!annotation$gene_id %in% tmp$gene,] + offset

    close2start <- tmp %>% dplyr::filter(start_range - offset < 0)
    if (nrow(close2start) != 0) {
        start(annotation[annotation$gene_id %in% close2start$gene,]) <- 1
        end(annotation[annotation$gene_id %in% close2start$gene,]) <-
           end(annotation[annotation$gene_id %in% close2start$gene,]) + offset
    }

    close2end <- tmp %>% dplyr::filter(end_range + offset > length)
    if (nrow(close2end) != 0) {
        end(annotation[annotation$gene_id %in% close2end$gene,]) <-
            close2end$length
        start(annotation[annotation$gene_id %in% close2end$gene,]) <-
            start(annotation[annotation$gene_id %in% close2end$gene,]) - offset
    }
    annotation
}

# define mouse chromosome ids
chr <- paste("chr", c(1:19, "M", "X", "Y"), sep="")

## Load annotations ####
hub <- AnnotationHub()
# |Organism: Mus musculus
# |taxonomy_id: 10090
# |genome_build: GRCm38
# |DBSCHEMAVERSION: 2.1
# | No. of genes: 56289.
# | No. of anno_transcriptss: 144726.
anno <- hub[["AH78811"]]

## Change seqlevels style to UCSC ie chr1 insted of 1; prevent warning by setting
## option for missing names
options(ensembldb.seqnameNotFound = "ORIGINAL")
seqlevelsStyle(anno) <- "UCSC"

# get genes, transcripts and exons
anno_genes <- genes(anno, filter=SeqNameFilter(chr))
anno_transcripts <- transcripts(anno, filter=SeqNameFilter(chr))
anno_exons <- exonsBy(anno, by='gene', filter=SeqNameFilter(chr))

# get genome data
genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10

# read TSS data
m5pseq <- data.table::fread(args$mouse, sep=",", header = TRUE,
                            stringsAsFactors=FALSE, data.table=FALSE)
fantom <- data.table::fread(args$fantom, sep=",", header = TRUE,
                            stringsAsFactors=FALSE, data.table=FALSE)

################
## analysis ####
################

## Generate feature table for annotation ####
features <- generate_features(anno_exons, anno_genes, anno_transcripts)
saveRDS(features, file.path(args$odir, "features_mouse_mm10.rds"))

## Pre-process tss data ####
total_m5pseq <- preprocess_feature_table(m5pseq, fantom, anno_genes)
total_m5pseq <- dplyr::rename(total_m5pseq, start=Position, Fantom5=response,
                       BarcodeCount=BarcodeCount.x)

total_fantom <- preprocess_feature_table(fantom, m5pseq, anno_genes)
total_fantom <- dplyr::rename(total_fantom, start=Position, m5pseq=response,
                       BarcodeCount=BarcodeCount.x)

## Annotate tss data with features ####
features_fantom <- collect_features(total_fantom, features, genome)
features_m5pseq <- collect_features(total_m5pseq, features, genome)

write.table(features_m5pseq,
            file=file.path(args$odir, "features_5pseq.csv"),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(features_fantom,
            file=file.path(args$odir, "features_fantom5.csv"),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
