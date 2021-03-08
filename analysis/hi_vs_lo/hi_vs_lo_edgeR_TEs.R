library(edgeR)

#################################################################
# Differential expression with edgeR
#################################################################

## Import count matrix

data = read.table("hi_vs_lo.cntTable",header=T,row.names=1)

## Subset into TE count matrix 'TE_data'

TE_data = data[grepl("^(?!ENSG).*$",rownames(data), perl = TRUE),]

## Create DGEList with filtering

group = c('hi', 'hi', 'hi', 'lo', 'lo', 'lo')
patient = c(214, 221, 226, 214, 221, 226)
y = DGEList(counts = TE_data, group = group)

y$samples$patient = patient

keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=FALSE]

y = calcNormFactors(y)

design = model.matrix(~patient +group)
y = estimateDisp(y,design)

fit <- glmQLFit(y, design)
tr <- glmTreat(fit, lfc=1)
topTags(tr, n = 20, p.value = 0.1)

#################################################################
#################################################################
# Data visualization
#################################################################
#################################################################