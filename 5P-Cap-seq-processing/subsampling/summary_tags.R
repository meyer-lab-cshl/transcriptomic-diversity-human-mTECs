log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#################
## libraries ####
#################

library(tidyverse)

################
## analysis ####
################


files <- list.files(snakemake@wildcards[['dir']],
                    pattern="tagCount.csv", recursive=TRUE,
                    include.dirs=TRUE, full.names=TRUE)

tags <- sapply(files, function(x) {
                   tmp <- read_csv(x)
                   tmp$reads <- sub(".*reads(.*)/ctss.*", "\\1", x)
                   return(tmp)
    }) %>%
    bind_rows()

tag_counts <- tags %>%
    group_by(reads) %>%
    select(starts_with("pt")) %>%
    summarise_all(sum)

write_csv(tag_counts, snakemake@output[[1]])

sessionInfo()
