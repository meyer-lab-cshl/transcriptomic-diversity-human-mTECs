snakemake -s processing_alignment_human.smk \
    --dag ~/data/tss/human/2_alignments/pt221-lo_5prime_nz_readdepth.combined.ucsc.bedgraph \
    ~/data/tss/combined/3_tss_data/raw_positions/all_mTECs_fwd.positions.csv |  dot -Tpdf > processing_alignment_human_dag.pdf
