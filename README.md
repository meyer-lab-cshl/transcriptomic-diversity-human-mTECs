# TE_thymus

TE annotation downloaded from http://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/

## Workflow:

1. Align reads with STAR
2. Quantify gene and locus-level TE expression with TElocal
3. Differential expression analysis with DESEq2
    + Prefiltering: minimum of one read per gene/TE
    + FDR: 0.1
    + Independent filtering: off