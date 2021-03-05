library(ggplot2)

##This is a set of plot templates for visualizing DESeq2 output. It assumes
##the following assignments:
##
## DESeq dataset object: dds
## Differential expression results: res
## Transformed counts: vs_dds

## PCA

plotPCA(vs_dds, intgroup = 'group') + theme_bw()

## Count plots for individual genes

par(mfrow=c(2,2))
plotCounts(dds, gene = 'ENSG00000160224.17', intgroup = 'group') ## AIRE
plotCounts(dds, gene = 'ENSG00000121594.12', intgroup = 'group') ## CD80
plotCounts(dds, gene = 'ENSG00000153266.13', intgroup = 'group') ## FEZF2
plotCounts(dds, gene = '', intgroup = 'group') ## FEZF2

## Volcano

df2 = filter(df, Class == 'protein_coding')

library(EnhancedVolcano)

EnhancedVolcano(df2,
                lab = df2$GeneSymbol,
                x = 'log2FoldChange',
                selectLab = c('AIRE', 'CD80', 'FEZF2'),
                y = 'padj',
                title = 'mTEC-lo vs mTEC-hi gene expression',
                pCutoff = 0.1,
                drawConnectors = TRUE)