#################
## libraries ####
#################

library(tidyverse)
library(plyr)
library(optparse)
library(data.table)
library(ggpointdensity)
library(viridis)

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
    args$fantomdir <- "~/data/tss/mouse/fantom/tss/isoforms"
    args$mousedir <- "~/data/tss/mouse/5Pseq/tss/isoforms"
    args$odir <- "~/data/tss/combined/tss/isoforms"
    args$suffix <- ".isoforms.csv"
}

fantomfiles <- list.files(path=args$fantomdir, pattern=args$suffix,
                          full.names=TRUE)

mousefiles <- list.files(path=args$mousedir, pattern=args$suffix,
                          full.names=TRUE)

# Iterate over all counts files in directory
fantomlist <- lapply(fantomfiles, function(s) {
    dat <- read_csv(s)
})

mouselist <- lapply(mousefiles, function(s) {
    dat <- read_csv(s)
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

dat_fantom_agg <- dat_fantom %>%
    group_by(geneID, Position) %>%
    summarise(BarcodeCount=mean(BarcodeCount)) %>%
    ungroup

# Combine counts from all replicates of 5Pseq mESC 46C samples
dat_mouse_agg <- dat_mouse %>%
    group_by(geneID, Position) %>%
    summarise(BarcodeCount=mean(BarcodeCount)) %>%
    ungroup


write_csv(dat_mouse_agg,
            file=file.path(args$odir, "all_mouse_mESC.isoforms.csv"))
write_csv(dat_fantom_agg,
          file.path(args$odir, "all_mouse_ES_46C_fantoms.isoforms.csv"))

# Combined 5Pseq and fantom counts
total <- full_join(dat_mouse_agg, dat_fantom_agg, by=c("geneID", "Position"))
total <- dplyr::rename(total,
                       c("CountInternal"="BarcodeCount.x",
                         "CountFantom"="BarcodeCount.y"))
write_csv(total, file.path(args$odir, "combined_mouse_ESC_5Pseq_fantom.csv"))

summary_total <- total %>%
    summarise(fantom_only=sum(is.na(CountInternal)),
              internal_only=sum(is.na(CountFantom)),
              both=sum(!(is.na(CountInternal) | is.na(CountFantom))))
write_csv(summary_total, file.path(args$odir, "comparison_tss_5Pseq_fantom.csv"))


# For visualisation on log scale:
total[is.na(total)] <- 1

total <- dplyr::filter(total, !is.na(CountInternal), !is.na(CountFantom))
tmp <- total %>%
    mutate(CountInternal=CountInternal/sum(total$CountInternal),
           CountFantom=CountFantom/sum(total$CountFantom))

tmp <- dplyr::filter(tmp, CountInternal > 0.00001,
                     CountFantom > 0.00001)

p <- ggplot(total, aes(x= CountInternal, y= CountFantom)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="Count mESC 5Pseq data",
       y= "Count mESC Fantom5 data") +
  scale_x_log10() +
  scale_y_log10() +
  coord_fixed() +
  theme_bw()
ggsave(plot=p, file.path(args$odir, "plots",
                         "counts_fantom_vs_5Pseq_per_position_density.pdf"))

