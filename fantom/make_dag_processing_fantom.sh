snakemake -s processing_fantom_mouse.smk \
    --rulegraph ~/data/tss/mouse/fantom/tss/combined/all_mESC_46C.positions.csv \
    ~/data/tss/mouse/fantom/tss/isoforms/ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1.isoforms.csv \
    |  dot -Tpng > dag/processing_fantom_dag.png

snakemake -s processing_fantom_human.smk \
    --rulegraph ~/data/tss/human/fantom/tss/combined/all_tissues.positions.csv \
    ~/data/tss/human/fantom/tss/isoforms/esophagus_adult_pool1.CNhs10620.10015-101C6.isoforms.csv \
    |  dot -Tpng > dag/processing_fantom_dag.png

#~/data/tss/combined/fantom_liftover_human_mouse.pdf 
