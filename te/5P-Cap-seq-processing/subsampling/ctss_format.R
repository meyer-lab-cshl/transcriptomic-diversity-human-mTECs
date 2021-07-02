#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

#################
## libraries ####
#################

library(limma)

################
## analysis ####
################

df <- read.csv(snakemake@input[['tpm']])
df_tagCounts <- read.csv(snakemake@input[['tagcounts']])


# Limma batch correction
# Define column order, must match batch order
samples_low <- c('pt212_lo','pt221_lo','pt226_lo','pt87_lo','pt214_lo')
samples_high <- c('pt212_hi','pt221_hi','pt226_hi','pt87_hi','pt214_hi')
samples <- c(samples_low, samples_high)
batch <- c(1, 2, 2, 1, 1, 1, 2, 2, 1, 1)
names(batch) <- samples

current_samples <- colnames(df_tagCounts)[grepl("pt", colnames(df_tagCounts))]
batch <- batch[names(batch) %in% current_samples]

if (all(c(1,2) %in% batch)) {
    df_batch <- removeBatchEffect(log2(df[,current_samples] + 0.01), batch)
    df_batch <- 2^df_batch
} else {
    df_batch <- df
}

# Prepare output for Paraclu
features <- c('chr', 'strand', 'pos')
for (subject in current_samples) {
    x <- df[,features]
    x$TPM <- df_batch[,subject]
    x$tagCounts <- df_tagCounts[,subject]
    subject <- gsub("_", "-", subject)
    write.table(x[x$tagCounts>0,],
                paste(snakemake@wildcards[['dir']], "/subsampling/reads",
                      snakemake@wildcards[['reads']], "/CAGEr_out/DFs/",
                      subject, "_Counts_and_TPM.txt", sep=""),
                row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    write.table(x[x$tagCounts>0, c(features,"TPM")],
                paste(snakemake@wildcards[['dir']], "/subsampling/reads",
                      snakemake@wildcards[['reads']], "/for_Paraclu/",
                      subject, "_TPM_for_Paraclu.txt", sep=""),
            row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}
