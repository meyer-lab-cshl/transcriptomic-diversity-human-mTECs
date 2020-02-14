snakemake -s processing_alignment.smk \
    --dag ~/data/tss/mouse/5Pseq/bedgraphs/mESC1_1_5prime_nz_readdepth.combined.bedgraph \
    ~/data/tss/mouse/5Pseq/tss/combined/all_mESCs.positions.csv |  dot -Tpng > dag/processing_alignment_dag.png
