#!/usr/bin/env Rscript

library(BSgenome.Hsapiens.UCSC.hg38)
library(CAGEr)
library(limma)

# Input Bam Files
analysis_dir <- "/grid/meyer/home/jacarter/TSS/CTSS/"
bam_files <- list.files("/grid/meyer/home/jacarter/TSS/5PSeq_Data",
                        full.names=TRUE, pattern="*.bam")
i=0
for (val in basename(bam_files)) {
    x <-strsplit(val, "_")
    if (i==0) {
        file_names <- x[[1]][1]
    } else {
        file_names <- c(file_names,x[[1]][1])
    }
    i=i+1
}

# Create CAGEr exp
ce <- CAGEexp(genomeName     = "BSgenome.Hsapiens.UCSC.hg38",
              inputFiles     = bam_files,
              inputFilesType = "bamPairedEnd",
              sampleLabels   = sub("-", "_", file_names)
)

# Call CTSS
getCTSS(ce, removeFirstG = FALSE,  correctSystematicG = FALSE)

pdf(file=paste(analysis_dir, 'PowerLaw.pdf', sep=''), width=7, height=7)
plotReverseCumulatives(ce, fitInRange = c(5, 2000), onePlot = TRUE)
dev.off()

# Power Law normalization -> tags per million
normalizeTagCount(ce, method = "powerLaw",  alpha = 1.15, T = 1*10^6)

pdf(file=paste(analysis_dir, 'PowerLaw_Normalized.pdf', sep=''),
    width=7, height=7)
plotReverseCumulatives(ce, values='normalized', fitInRange = c(5, 2000),
                       onePlot = TRUE)
dev.off()

# Save CTSS across samples with TPM and raw tag counts
df <- CTSSnormalizedTpm(ce)
write.csv(df, paste(analysis_dir, "CAGEr_out/TPM.csv", sep=''))

df_tagCounts <- CTSStagCount(ce)
write.csv(df_tagCounts, paste(analysis_dir, "CAGEr_out/tagCount.csv", sep=''))

# Limma batch correction
# Define column order, must match batch order
samples_low <- c('pt212_lo','pt221_lo','pt226_lo','pt87_lo','pt214_lo')
samples_high <- c('pt212_hi','pt221_hi','pt226_hi','pt87_hi','pt214_hi')
samples <- c(samples_low,samples_high)

batch <- c(1, 2, 2, 1, 1, 1, 2, 2, 1, 1)
y2 <- removeBatchEffect(log2(df[,samples]+0.01), batch)
y2 <- 2^y2

# Prepare output for Paraclu
features=c('chr','strand','pos')
for (subject in samples) {
    x <- df[,features]
    x$TPM <- y2[,subject]
    x$tagCounts <- df_tagCounts[,subject]
    write.table(x[x$tagCounts>0,],
                paste(analysis_dir, "CAGEr_out/DFs/", subject,
                      "_Counts_and_TPM", sep=""),
                row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    write.table(x[x$tagCounts>0, c(features,"TPM")],
                paste(analysis_dir, "for_Paraclu/", subject,
                      "_TPM_for_Paraclu.txt", sep=""),
            row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}
