log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#################
## libraries ####
#################

library(BSgenome.Hsapiens.UCSC.hg38)
library(CAGEr)
library(limma)

sessionInfo()
################
## analysis ####
################

cat(unlist(snakemake@input))
# Create CAGEr exp
ce <- CAGEexp(genomeName     = "BSgenome.Hsapiens.UCSC.hg38",
              inputFiles     = unlist(snakemake@input),
              inputFilesType = "bamPairedEnd",
              sampleLabels   = sub("-", "_", snakemake@params[['samples']])
)

# Call CTSS
getCTSS(ce, removeFirstG = FALSE,  correctSystematicG = FALSE)

pdf(file=snakemake@output[['powerlaw']], width=7, height=7)
plotReverseCumulatives(ce, fitInRange = c(5, 2000), onePlot = TRUE)
dev.off()

# Power Law normalization -> tags per million
normalizeTagCount(ce, method = "powerLaw",  alpha = 1.15, T = 1*10^6)

pdf(file=snakemake@output[['powerlawnorm']], width=7, height=7)
plotReverseCumulatives(ce, values='normalized', fitInRange = c(5, 2000),
                       onePlot = TRUE)
dev.off()

# Save CTSS across samples with TPM and raw tag counts
df <- CTSSnormalizedTpm(ce)
write.csv(df, snakemake@output[['tpm']])

df_tagCounts <- CTSStagCount(ce)
write.csv(df_tagCounts, snakemake@output[['tagcounts']])

