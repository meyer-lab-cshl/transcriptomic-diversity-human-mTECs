log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

#reads <- c(0.5, 1:10, 15, 20, 25, 30) * 10^5
reads <- c(seq(3,39, 3), 45, 60, 75, 90, 105, 120) * 10^6

overview_reads <- read_delim(snakemake@input[[1]],
                             delim="\t", col_names = c("id", "total_reads")) %>%
    mutate(id =(str_squish(id)),
           total_reads = as.numeric(str_squish(total_reads)))

subsampling <- tibble(id=rep(overview_reads$id, each=length(reads)),
                      subsample_reads=rep(reads, nrow(overview_reads))) %>%
    left_join(overview_reads, by="id") %>%
    mutate(proportion = subsample_reads/total_reads) %>%
    filter(proportion < 1)

write_delim(subsampling, snakemake@output[[1]], delim="\t")
