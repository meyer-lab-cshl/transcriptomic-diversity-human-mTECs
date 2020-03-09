#################
## libraries ####
#################
library(tidyverse)
library(optparse)

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--dir"), action="store", dest="dir",
                type="character", help="Path to output directory
                data [default: %default].", default=NULL),
    make_option(c("--human"), action="store", dest="human",
                type="character", help="Path to file with human fantom liftover
                summary [default: %default].", default=NULL),
    make_option(c("--mouse"), action="store", dest="mouse",
                type="character", help="Path to file with mouse fantom liftover
                summary [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$human <- "~/data/tss/human/fantom/bed/GRCh38/all_tissues.count.unmapped.txt"
    args$mouse <- "~/data/tss/mouse/fantom/bed/GRCm38/all_mESC_C46.count.unmapped.txt"
    args$dir <- "~/data/tss/combined"

}

#############
## data  ####
#############

human <- read_delim(args$human, col_names=c("sample", "mapped", "unmapped"),
                    delim="\t")
mouse <- read_delim(args$mouse, col_names=c("sample", "mapped", "unmapped"),
                    delim="\t")

################
## analysis ####
################

## combine datasets ####
human <- human %>%
    mutate(species="human")

mouse <- mouse %>%
    mutate(species="mouse")

combined <- bind_rows(human, mouse) %>%
    mutate(ppt=unmapped/mapped)
write_csv(combined, file.path(args$dir, "fantom_liftover_human_mouse.csv"))

## display overview ####
p <- ggplot(combined, aes(y=ppt*100, x=species, color=species))
p <- p + geom_boxplot() +
    scale_color_brewer(type="qual", palette="Set1", guide=FALSE) +
    labs(y="Fantom5 liftover unmapped [%]",
         x="Data set") +
    theme_bw()

ggsave(plot=p, file.path(args$dir, "fantom_liftover_human_mouse.pdf",
    height=7, width=7, unit="cm")
