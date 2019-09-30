#################
## Libraries ####
#################

library(data.table)
library(optparse)


#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--directory"), action="store", dest="directory",
               type="character", help="Path to input directory
               [default: %default].", default=NULL),
    make_option(c("-s", "--suffix"), action="store", dest="suffix",
               type="character", help="common suffix of files to read
               [default: %default].", default=NULL),
    make_option(c("-o", "--ofile"), action="store", dest="ofile",
               type="character", help="Path to output file [default: %default].",
               default=NULL),
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
    args$directory <- "~/data/common/tss/3_tss_data/raw_positions/single_samples"
    args$ofile <- "~/data/common/tss/3_tss_data/raw_positions/all_mTECs.positions.csv"
    args$suffix <- "_fwd.positions.csv"
    args$verbose <- TRUE
}

samplefiles <- list.files(path=args$directory, pattern=args$suffix,
                          full.names=TRUE)

# Iterate over all counts files in directory
datalist <- lapply(samplefiles, function(s) {
    # Read counts per sample
    dat <- fread(s, header=TRUE, stringsAsFactors=FALSE)
})

dat_all <- do.call(rbind, datalist)

################
## Analysis ####
################

if (args$verbose) message("aggregating")
dat_agg_all <- dat_all[,lapply(.SD, sum),
                       by=.(chrom, geneID, strand, start, end),
                       .SDcols=c("BarcodeCount")]

if (args$verbose) message("writing out")
write.table(setDF(dat_agg_all), file=args$ofile, row.names=FALSE, na="",
            col.names=TRUE, quote=FALSE, sep=",")
