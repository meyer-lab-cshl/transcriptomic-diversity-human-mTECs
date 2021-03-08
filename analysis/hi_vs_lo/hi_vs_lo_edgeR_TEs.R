library(edgeR)

#################################################################
# Differential expression with edgeR
#################################################################

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)

## Subset into TE count matrix 'TE_data'

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Create DGEList with filtering

group <- c(1, 1, 1, 2, 2, 2)
y = DGEList(counts = TE_data, group = group)

keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=FALSE]

y = calcNormFactors(y)

design = model.matrix(~group)
y = estimateDisp(y,design)

et <- exactTest(y)
topTags(et)