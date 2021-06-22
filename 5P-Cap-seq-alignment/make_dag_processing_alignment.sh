snakemake -s processing_alignment.smk \
    --rulegraph  ~/data/tss/human/5Pseq/bedgraphs/pt221-lo_5prime_nz_readdepth.combined.ucsc.bedgraph \
    |  dot -Tpng > dag/processing_alignment_dag.png
