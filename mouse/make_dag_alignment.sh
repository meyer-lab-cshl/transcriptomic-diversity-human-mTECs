snakemake -s alignment.smk \
    --rulegraph ~/data/tss/mouse/5Pseq/alignments/mESC1_1_Aligned.sortedByCoord.out.bam \
    ~/data/tss/mouse/5Pseq/multiqc/multiqc_report.html \
    ~/data/tss/mouse/5Pseq/alignment_qc/multiqc_report.html | dot -Tpng > dag/alignment_dag.png
