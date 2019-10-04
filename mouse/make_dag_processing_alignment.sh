snakemake -s processing_alignment_mouse.smk \
    --dag ~/data/tss/mouse/2_alignments/mESC1_1_5prime_nz_readdepth.combined.bedgraph \
    ~/data/tss/combined/3_tss_data/raw_positions/all_mESCs.positions.csv |  dot -Tpdf > processing_alignment_mouse_dag.pdf
