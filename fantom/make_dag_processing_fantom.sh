snakemake -s processing_fantom.smk \
    --dag ~/data/tss/combined/3_tss_data/raw_positions/all_tissues.hg19.ctss.positions.csv \
         ~/data/tss/combined/3_tss_data/raw_positions/all_mESC_46C.mm9.ctss.positions.csv |  dot -Tpdf > processing_fantom_dag.pdf
