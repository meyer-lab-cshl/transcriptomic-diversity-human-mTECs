snakemake -s processing_alignment_human.smk \
    --dag /mnt/grid/meyer/hpc/home/data/common/tss/2_alignments/pt221-lo_5prime_nz_readdepth.combined.ucsc.bedgraph \
    /mnt/grid/meyer/hpc/home/data/common/tss/3_tss_data/raw_positions/single_samples/pt221-lo_fwd.positions.csv |  dot -Tpdf > processing_alignment_human_dag.pdf
