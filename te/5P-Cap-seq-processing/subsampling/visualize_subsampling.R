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

textsize <- 6
titlesize <- 7

tagcounts <- read_csv(snakemake@input[['tags']]) %>%
    pivot_longer(-reads, names_to="sample", values_to="tags") %>%
    separate(sample, into=c("id", "type"), sep="_", remove=FALSE) %>%
    drop_na

p <- ggplot(tagcounts,
            aes(x=reads/10^6, y=tags/10^6, group=sample, color=id, shape=type)) +
    geom_point() +
    geom_line() +
    scale_color_brewer(type="qual", palette = "Dark2") +
    labs(x="Number of molecules [x10^6]",
         y="Number of TSS [x10^6]",
         shape = "Stage",
         color = "Sample ID") +
    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    cowplot::theme_cowplot() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box="vertical",
          axis.text = element_text(size=textsize),
          legend.text = element_text(size=textsize),
          axis.title = element_text(size=titlesize),
          legend.title = element_text(size=titlesize))
ggsave(plot=p, filename = snakemake@output[['tags']],
       height = 9, width=6,  units = "cm")


clustercounts <- read_csv(snakemake@input[['clusters']]) %>%
    separate(sample, into=c("id", "type"), sep="-", remove=FALSE)


p <- ggplot(clustercounts,
       aes(x=reads/10^6, y=clusters/10^4, group=sample, color=id, shape=type)) +
    geom_point() +
    geom_line() +
    scale_color_brewer(type="qual", palette = "Dark2") +
    labs(x="Number of molecules [x10^6]",
         y="Number of TSR [x10^4]",
         shape = "Stage",
         color = "Sample ID") +
    guides(color=guide_legend(nrow=2, byrow=TRUE)) +
    cowplot::theme_cowplot() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box="vertical",
          axis.text = element_text(size=textsize),
          legend.text = element_text(size=textsize),
          axis.title = element_text(size=titlesize),
          legend.title = element_text(size=titlesize))
ggsave(plot=p, filename = snakemake@output[['clusters']],
       height = 9, width=6,  units = "cm")


