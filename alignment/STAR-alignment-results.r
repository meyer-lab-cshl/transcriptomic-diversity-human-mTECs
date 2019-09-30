#################
## Libraries ####
#################

library(data.table)
library(tidyverse)
library(optparse)

parse_star <- function(starfile, suffix) {
        sample <- gsub(paste(".*/(pt.*)", suffix, sep=""), "\\1", starfile)
        df <- data.table::fread(starfile, sep="\t", data.table=FALSE, fill=TRUE,
        stringsAsFactors=FALSE, skip=3)

        total_reads <- as.numeric(df[grepl("^Number of input", df[,1]),2])
        uniquely_mapped <- as.numeric(df[grepl("^Uniquely mapped reads number", df[,1]),2])
        multiple_mapped <- as.numeric(df[grepl("^Number of reads mapped to multiple", df[,1]),2])
        many_mapped <- as.numeric(df[grepl("^Number of reads mapped to too many", df[,1]),2])
        short_unmapped <- as.numeric(df[grepl("^Number of reads unmapped: too short", df[,1]),2])
        mismatches_unmapped <- as.numeric(df[grepl("^Number of reads unmapped: too many mismatches", df[,1]),2])
        other_unmapped <- as.numeric(df[grepl("^Number of reads unmapped: other", df[,1]),2])

        all_multiple_mapped <- multiple_mapped + many_mapped
        all_unmapped <- short_unmapped + mismatches_unmapped + other_unmapped

        reads <- data.frame(patient=sample, total_reads=total_reads,
                            uniquely_mapped=uniquely_mapped,
                            multi_mapped=all_multiple_mapped,
                            unmapped=all_unmapped, stringsAsFactors=FALSE)
        return(reads)
}

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--directory"), action="store", dest="directory",
               type="character", help="Path to directory with alignent files
               [default: %default].", default=NULL),
    make_option(c("-s", "--suffix"), action="store", dest="suffix",
               type="character", help="common suffix of files to read
               [default: %default].", default=NULL),
    make_option(c("-v", "--verbose"), action="store_true", dest="verbose",
               type="logical", help="Print progress to stdout
               [default: %default].", default=FALSE),
    make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$directory <- "~/data/common/tss/2_alignments"
    args$suffix <- "_Log.final.out"
    args$verbose <- TRUE
}

samplefiles <- list.files(path=args$directory, pattern=args$suffix,
                          full.names=TRUE)

################
## Analysis ####
################

# Iterate over all counts files in directory
datalist <- lapply(samplefiles, function(s) {
    all_reads <- parse_star(s, args$suffix)
})

reads <- do.call(rbind, datalist)

write.table(reads, file.path(args$directory, "STAR_summary.csv"),
        col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")

reads <- reads %>%
            pivot_longer(c(-patient, -total_reads), names_to="type",
                         values_to="count") %>%
            mutate(percent=count*100/total_reads) %>%
            separate(patient, into=c("sample", "mtecs"), sep="-")
reads$sample <- factor(reads$sample,
                       levels=paste("pt", c('87','212','214','221','226'),
                                    sep=""))

reads$type <- factor(reads$type,
                       levels=c("uniquely_mapped", "multi_mapped", "unmapped"),
                       labels=c("uniquely mapped", "multi mapped", "unmapped"))

p <- ggplot(reads, aes(x=sample, y=percent, fill=type))
p <- p + geom_bar(stat='identity', position='dodge') +
         labs(x ="Patient sample", y ="Percentage of reads [%]" ) +
         scale_fill_brewer(type="qual", palette=2) +
         guides(fill=guide_legend(title="STAR category")) +
         theme_bw()


ggsave(plot=p, file="STAR_summary.pdf", path = args$directory)
