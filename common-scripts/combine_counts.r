#################
## Libraries ####
#################

library(data.table)
library(optparse)


#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-i", "--indir"), action="store", dest="indir",
               type="character", help="Path to input directory
               [default: %default].", default=NULL),
    make_option(c("-o", "--odir"), action="store", dest="odir",
               type="character", help="Path to output directory
               [default: %default].", default=NULL),
    make_option(c("-s", "--suffix"), action="store", dest="suffix",
               type="character", help="common suffix of files to read
               [default: %default].", default=NULL),
    make_option(c("-f", "--fantom"), action="store_true", dest="verbose",
               type="logical", help="Fantom data; process thymus separately
               [default: %default].", default=FALSE),
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
    args$indir <- "~/data/human/tss/3_tss_data/raw_positions"
    args$odir <- "~/data/common/tss/3_tss_data/raw_positions"
    args$suffix <- "_fwd.positions.csv"
    args$fantom <- TRUE
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

if(args$fantom) {
    wo_thymus <- datalist[!grepl("thymus", datalist)]
    dat_wo_thymus <- do.call(rbind, wo_thymus)
}


################
## Analysis ####
################

if (args$verbose) message("aggregating")
dat_agg_all <- dat_all[,lapply(.SD, sum),
                       by=.(chrom, geneID, strand, start, end),
                       .SDcols=c("BarcodeCount")]

if (args$verbose) message("writing out")
write.table(setDF(dat_agg_all),
            file=file.path(args$odir, paste("all_tissues", args$suffix, sep="")),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

if(args$fantom) {
    if (args$verbose) message("aggregating w/o thymus samples")
    dat_agg_wo_thymus <- dat_wo_thymus[,lapply(.SD, sum),
                                       by=.(chrom, geneID, strand, start, end),
                                       .SDcols=c("BarcodeCount")]

    if (args$verbose) message("writing out aggregate w/o thymus")
    write.table(setDF(dat_agg_all),
                file=file.path(args$odir, paste("all_tissues_wo_thymus",
                                                args$suffix, sep="")),
                row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
}
