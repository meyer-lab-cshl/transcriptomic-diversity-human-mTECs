snakemake -s tss.smk \
    --rulegraph ~/data/tss/human/5Pseq/tss/summary/all_samples_consensus_sites.csv \
    ~/data/tss/human/5Pseq/tss/isoforms/pt212-hi.isoforms.csv |  dot -Tpng > dag/tss_dag.png
