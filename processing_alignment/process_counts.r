#################
## Libraries ####
#################

library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(optparse)


#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-o", "--odir"), action="store", dest="outdir",
               type="character", help="Path to output directory
               [default: %default].", default=NULL),
    make_option(c("-s", "--sample"), action="store", dest="sample",
               type="character", help="Name of sample [default: %default].",
               default=NULL),
    make_option(c("-i", "--ifile"), action="store", dest="ifile",
               type="character", help="Path to input file [default: %default].",
               default=NULL),
    make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$odir <- "~/data/common/tss/3_tss_data"
    args$sample <- "pt214-lo"
    args$ifile <- paste("~/data/common/tss/2_alignments/collapsing_fwd/",
                        "pt221-lo_collapsed_countsPerMol.tsv", sep="")
}

## Define chromosome IDs
chr <- paste ("chr", c(1:22, "M", "X", "Y"), sep="")

## Load gene object and components
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb) <- chr

# get gene ids, expand range by 100, order by gene_id
tx_genes <- genes(txdb)
tx_genes <- tx_genes + 100
tx_genes <- tx_genes[order(as.numeric(tx_genes$gene_id))]

# Read counts per molecule
dat <- data.table::fread(args$ifile, header=TRUE, data.table=FALSE,
                         stringsAsFactors=FALSE)

################
## Analysis ####
################

# drop non-standard chromosomes and IVT
dat <- dat[grepl("^chr", dat$chr),]

# drop molbc and counts
dat <- dplyr::select(dat, -molbc, -count)
dat$poscount <-  1

# calculate real count based on positions
dat_agg <- aggregate(dat$poscount, by=list(dat$chr, dat$pos, dat$strand),
                     FUN=sum)

# rename columns
dat_agg <- rename(dat_agg, c("Group.1"="chromosome", "Group.2"="start",
                             "Group.3"="strand", "x"="count"))
dat_agg$end <- as.integer(dat_agg$start + 1)
data_agg <- dplyr::select(dat_agg, chromosome, start, strand, end, count)

# normalise counts based on total read count
dat_agg$count <- dat_agg$count/(sum(dat_agg$count))*10000000

## Format into bedgraph ####
dat_bedgraph <- dat_agg
minus_strand <- dat_bedgraph$strand == "-"
dat_bedgraph$count[minus_strand] <- dat_bedgraph$count[minus_strand]*-1
dat_bedgraph <- dat_bedgraph[order(dat_bedgraph$chr, dat_bedgraph$start),]
dat_bedgraph <- data.frame("chrom"=dat_bedgraph$chromosome,
                           "start"=dat_bedgraph$start,
                           "end"=dat_bedgraph$end,
                           "dataValue"=dat_bedgraph$count)
dat_bedgraph <- dat_bedgraph[complete.cases(dat_bedgraph),]
write.table(dat_bedgraph,
            file=file.path(args$odir, "raw_positions", "bedgraphs",
                           paste(args$sample, ".bedgraph", sep="")),
            row.names=FALSE, na="", col.names=FALSE, quote=FALSE, sep="\t")

## find Overlaps with genes by comparing genomic ranges ####
dat_ranges <- makeGRangesFromDataFrame(dat_agg, keep.extra.columns = TRUE)
dat_agg$overlaps  <- findOverlaps(dat_ranges, tx_genes, select = 'first')

# lookup entrezId with overlap number
dat_agg$regions= with(dat_agg, names(tx_genes)[overlaps])

# set intergenic_$start/1000 instead of NA, to prevent loosing information
na.region <- dplyr::select(dat_agg, regions) %>% is.na
dat_agg$regions[na.region] <- paste("intergenic", dat_agg$chromosome[na.region],
                                    dat_agg$strand[na.region],
                                    as.integer(dat_agg$start[na.region]/1000),
                                    sep="_")
# Keep rows with no missing values only
dat_agg <- dat_agg[complete.cases(dat_agg),]

# format summary counts
dat_summary <- data.frame("geneID"=dat_agg$regions, "Position"=dat_agg$start,
                         "BarcodeCount"=dat_agg$count,
                         "ReadCount"=dat_agg$count, "PosFromAnno"=0, "Class"=0)
dat_summary <- dat_summary[with(dat_summary, order(geneID, Position)),]
write.table(dat_summary,
            file=file.path(args$odir, "raw_positions", "single_samples",
                           paste(args$sample, ".summary.counts",
                                            sep="")),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

# format summary positions
dat_positions <- data.frame("chrom"=dat_agg$chromosome, "geneID"=dat_agg$regions,
                         "strand"=dat_agg$strand, "start"=dat_agg$start,
                         "end"=dat_agg$end, "BarcodeCount"=dat_agg$count)

dat_positions <- dat_positions[with(dat_positions, order(geneID, start)),]
write.table(dat_positions,
            file=file.path(args$odir, "raw_positions", "single_samples",
                           paste(args$sample, ".positions.csv", sep="")),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")




