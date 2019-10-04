#################
## libraries ####
#################

library(plyr)
library(data.table)
library(ggplot2)
library(LSD)

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-m", "--mousedir"), action="store", dest="mousedir",
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
    args$fantomdir <- "~/data/tss/fantom/3_tss_data/raw_positions"
    args$mousedir <- "~/data/tss/mouse/3_tss_data/raw_positions"
    args$odir <- "~/data/tss/combined/3_tss_data/isoforms"
    args$suffix <- ".isoforms.csv"
}

fantomfiles <- list.files(path=args$fantomdir, pattern=args$suffix,
                          full.names=TRUE)

mousefiles <- list.files(path=args$mousedir, pattern=args$suffix,
                          full.names=TRUE)

# Iterate over all counts files in directory
fantomlist <- lapply(fantomfiles, function(s) {
    dat <- fread(s, header=TRUE, stringsAsFactors=FALSE)
})

mouselist <- lapply(mousefiles, function(s) {
    dat <- fread(s, header=TRUE, stringsAsFactors=FALSE)
})

dat_fantom <- do.call(rbind, fantomlist)
dat_mouse <- do.call(rbind, mouselist)

################
## analysis ####
################

# Combine counts from all replicates of fantom5 mESC 46C samples
dat_fantom_agg <- dat_fantom[, lapply(.SD, sum), by=.(geneID, Position),
                               .SDcols=c("BarcodeCount", "ReadCount",
                                         "PosFromAnno", "Class")]

# Combine counts from all replicates of 5Pseq mESC 46C samples
dat_mouse_agg <- dat_mouse[, lapply(.SD, sum), by=.(geneID, Position),
                               .SDcols=c("BarcodeCount", "ReadCount",
                                         "PosFromAnno", "Class")]

setDF(dat_fantom_agg)
setDF(dat_mouse_agg)

write.table(dat_mouse_agg,
            file=file.path(args$odir, "all_mouse_mESC.isoforms.csv"),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")
write.table(dat_fantom_agg,
            file=file.path(args$odir, "all_mouse_ES_46C_fantoms.isoforms.csv"),
            row.names=FALSE, na="", col.names=TRUE, quote=FALSE, sep=",")

# Combined 5Pseq and fantom counts
total <- merge(dat_mouse_agg, dat_fantom_agg, by=c("geneID", "Position"),
               all=TRUE)
total <- subset(total, select = -c(ReadCount.x, PosFromAnno.x, Class.x,
                                   ReadCount.y, PosFromAnno.y, Class.y))
total <- rename(total,c("BarcodeCount.x"="CountInternal",
                        "BarcodeCount.y"="CountFantom"))

pdf(file.path(args$odir, "plots", "Counts_fantom_vs_5Pseq_per_position.pdf")
heatscatter(total$CountInternal, total$CountFantom, log = "xy", 
            main="Counts per position",
            xlab="Count mESC 5Pseq data",
            ylab="Count mESC Fantom5 data")

ggplot(total, aes(x= CountInternal, y= CountFantom)) +
  geom_point() +
  labs(title="Counts per position",
       x="Count mESC 5Pseq data",
       y= "Count mESC Fantom5 data") +
  scale_x_log10() +
  scale_y_log10()
dev.off()

