# Endogenous retroviral elements in mTECs
Transposable element expression in medullary thymic epithelial cells.

## Workflow:

1. **Pre-processing with FastP**

2. **Align reads with STAR**
    + GENCODE GRCh38 gene annotation
    + sjdbOverhang: 100
    + winAnchorMultimapNMax: 200
    + outFilterMultimapNmax: 100
    
3. **Quantify gene and TE expression with TEcounts**
    + GENCODE GRCh38 Repeatmasker TE annotation (downloaded from http://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/)
    + GENCODE GRCh38 gene annotation
    
4. **Differential expression analysis with DESEq2**
    + Pre-filtering: minimum of 2 normalized reads per gene/TE
    + FDR: 0.1
    + Independent filtering: off