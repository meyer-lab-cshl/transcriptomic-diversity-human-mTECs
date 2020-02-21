snakemake -s processing_fantom.smk \
    --rulegraph ~/data/tss/mouse/fantom/tss/combined/all_mESC_46C.positions.csv \
    ~/data/tss/human/fantom/tss/combined/all_tissues.positions.csv |  dot -Tpng > dag/processing_fantom_dag.png
