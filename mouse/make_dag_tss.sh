snakemake -s tss.smk \
    --rulegraph ~/data/tss/mouse/5Pseq/tss/summary/all_samples_consensus_sites.csv \
    ~/data/tss/mouse/5Pseq/tss/isoforms/mESC1_1.isoforms.csv |  dot -Tpng > dag/tss_dag.png
