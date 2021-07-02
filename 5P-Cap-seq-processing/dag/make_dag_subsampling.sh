snakemake -s ../subsampling.smk \
    --rulegraph ~/data/tss/human/5Pseq/subsampling/subsampling_tags_summary.pdf | dot -Tpng > subsampling_dag.png
