#################
## Libraries ####
#################

library(tidyverse)
library(GenomicRanges)
library(optparse)
library(AnnotationHub)


extendGeneRange <- function(annotation, offset=100) {
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

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-o", "--odir"), action="store", dest="odir",
               type="character", help="Path to output directory
               [default: %default].", default=NULL),
    make_option(c("-s", "--sample"), action="store", dest="sample",
               type="character", help="Name of sample [default: %default].",
               default=NULL),
    make_option(c("-p", "--species"), action="store", dest="species",
               type="character", help="Name of species for gene mapping
               (mouse | human) [default: %default].",
               default=NULL),
    make_option(c("-n", "--normalise"), action="store_true", dest="normalise",
               type="character", help="Normalise counts by by total number of
               counts [default: %default].",
               default=FALSE),
    make_option(c("-i", "--ifile"), action="store", dest="ifile",
               type="character", help="Path to input file [default: %default].",
               default=NULL),
    make_option(c("-t", "--type"), action="store",
               dest="type", default=NULL, type="character",
               help="Specify data type (molbc|bed|bedgraph)
               [default: %default]."),
    make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$normalise <- TRUE
    args$type <- "bed"
    args$species <- "mouse"
    args$odir <- "~/data/tss/mouse/fantom/tss"
    args$sample <- "ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1"
    args$ifile <- paste("~/data/tss/mouse/fantom/bed/GRCm38/",
                        "ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1.bed", sep="")
}

## Load annotations ####
hub <- AnnotationHub()
if (args$species == "human") {
    # |Organism: Homo sapiens
    # |taxonomy_id: 9606
    # |genome_build: GRCh38
    # |DBSCHEMAVERSION: 2.1
    # | No. of genes: 67998.
    # | No. of transcripts: 250641.
    anno <- hub[["AH78783"]]
    chr <- paste ("chr", c(1:22, "M", "X", "Y"), sep="")
} else if(args$species == "mouse") {
    # |Organism: Mus musculus
    # |taxonomy_id: 10090
    # |genome_build: GRCm38
    # |DBSCHEMAVERSION: 2.1
    # | No. of genes: 56289.
    # | No. of transcripts: 144726.
    anno <- hub[["AH78811"]]
    chr <- paste ("chr", c(1:19, "M", "X", "Y"), sep="")
} else {
    stop("Species", args$species, "unknown")
}

## Change seqlevels style to UCSC ie chr1 insted of 1; prevent warning by setting
## option for missing names
options(ensembldb.seqnameNotFound = "ORIGINAL")
seqlevelsStyle(anno) <- "UCSC"

# get gene ids, expand range by 100
anno_genes <- genes(anno, filter=SeqNameFilter(chr))
anno_genes <- extendGeneRange(anno_genes, offset= 100)

################
## Analysis ####
################
if (args$type == "pos"){
    dat <- data.table::fread(args$ifile, header=TRUE, data.table=FALSE,
                             stringsAsFactors=FALSE)
    colnames(dat) <- c("chromosome", "strand", "start", "count")
} else if (args$type == "bedgraph") {
    dat <- data.table::fread(args$ifile, header=FALSE, data.table=FALSE,
                            stringsAsFactors=FALSE)
    colnames(dat) <- c('chromosome', 'start', 'end', 'count')
} else if (args$type == "bed") {
    dat <- data.table::fread(args$ifile, header=FALSE, data.table=FALSE,
                            stringsAsFactors=FALSE)
    colnames(dat) <- c("chromosome", "start", "end", "name", "count", "strand")
} else {
    stop("Data type", args$type, "unknown")
}

# drop non-standard chromosomes and IVT
dat_agg <- dat[grepl("^chr", dat$chromosome),]

if (args$type == 'pos') {
    dat_agg$end <- as.integer(dat_agg$start + 1)
    data_agg <- dplyr::select(dat_agg, chromosome, start, strand, end, count)
}

if (args$normalise) {
    # normalise counts based on total read count
    dat_agg$count <- dat_agg$count/(sum(abs(dat_agg$count)))*10000000
}

if (args$type == "pos" || args$type == "bed") {
    ## Format into bedgraph ####
    dat_bedgraph <- dat_agg
    minus_strand <- dat_bedgraph$strand == "-"
    dat_bedgraph$count[minus_strand] <- dat_bedgraph$count[minus_strand]*-1
    dat_bedgraph <- dat_bedgraph[order(dat_bedgraph$chromosome,
                                       dat_bedgraph$start),]
    dat_bedgraph <- data.frame("chrom"=dat_bedgraph$chromosome,
                               "start"=dat_bedgraph$start,
                               "end"=dat_bedgraph$end,
                               "dataValue"=dat_bedgraph$count)

    dat_bedgraph <- dat_bedgraph[complete.cases(dat_bedgraph),]
    write.table(dat_bedgraph,
                file=file.path(args$odir, "bedgraphs",
                               paste(args$sample, ".bedgraph", sep="")),
                row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
} else {
    dat_bedgraph <- dat_agg[complete.cases(dat_agg),]
    write.table(dat_agg,
                file=file.path(args$odir, "bedgraphs",
                               paste(args$sample, ".bedgraph", sep="")),
                row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
    dat_agg$strand <- "+"
    dat_agg$strand[sign(dat_agg$count) < 0] <-"-"
    dat_agg$count <- abs(dat_agg$count)
    dat_agg <- dplyr::select(dat_agg, chromosome, start, strand, end, count)
}

## find Overlaps with genes by comparing genomic ranges ####
dat_ranges <- makeGRangesFromDataFrame(dat_agg, keep.extra.columns = TRUE)
dat_agg$overlaps  <- findOverlaps(dat_ranges, anno_genes, select = 'first')

# lookup entrezId with overlap number
dat_agg$regions <- with(dat_agg, names(anno_genes)[overlaps])

# set intergenic_$start/1000 instead of NA, to prevent loosing information
na.region <- dplyr::select(dat_agg, regions) %>% is.na
dat_agg$regions[na.region] <- paste("intergenic", dat_agg$chromosome[na.region],
                                    dat_agg$strand[na.region],
                                    as.integer(dat_agg$start[na.region]/1000),
                                    sep="_")
# Keep rows with no missing values only
dat_agg <- dat_agg[complete.cases(dplyr::select(dat_agg, -overlaps)),]

# format summary counts
dat_summary <- data.frame(geneID=dat_agg$regions, Position=dat_agg$start,
                          BarcodeCount=dat_agg$count,
                          ReadCount=dat_agg$count, PosFromAnno=0, Class=0,
                          stringsAsFactors=FALSE)
dat_summary <- dat_summary[with(dat_summary, order(geneID, Position)),]
write.table(dat_summary,
            file=file.path(args$odir, "summary",
                           paste(args$sample, ".summary.counts.csv",
                                            sep="")),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

# format summary positions
dat_positions <- data.frame(chrom=dat_agg$chromosome, geneID=dat_agg$regions,
                            strand=dat_agg$strand, start=dat_agg$start,
                            end=dat_agg$end, BarcodeCount=dat_agg$count,
                            stringsAsFactors=FALSE)

dat_positions <- dat_positions[with(dat_positions, order(geneID, start)),]
write.table(dat_positions,
            file=file.path(args$odir, "raw_positions",
                           paste(args$sample, ".positions.csv", sep="")),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

# format summary ctts
dat_ctts <- data.frame(chrom=dat_agg$chromosome,
                            start=dat_agg$start,
                            strand=dat_agg$strand,
                            BarcodeCount=dat_agg$count,
                            stringsAsFactors=FALSE)

dat_ctts <- dat_ctts[with(dat_ctts, order(chrom, start)),]
write.table(dat_ctts,
            file=file.path(args$odir, "ctts",
                           paste(args$sample, ".ctts.txt", sep="")),
            row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")

# format summary bed
dat_bed <- data.frame(chrom=dat_agg$chromosome,
                      start=dat_agg$start,
                      end=dat_agg$end,
                      name=paste(dat_agg$chromosome, ":", dat_agg$start, "..",
                                 dat_agg$end, ",", dat_agg$strand, sep=""),
                      BarcodeCount=dat_agg$count,
                      strand=dat_agg$strand,
                      stringsAsFactors=FALSE)

dat_bed <- dat_bed[with(dat_bed, order(chrom, start)),]
write.table(dat_bed,
            file=file.path(args$odir, "ctts",
                           paste(args$sample, ".bed", sep="")),
            row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")
