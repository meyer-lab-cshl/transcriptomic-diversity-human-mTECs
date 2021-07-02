snakemake -s alignment.smk \
    --rulegraph ~/data/tss/human/5Pseq/fastq/pt221-lo_1.fastq \
    ~/data/tss/human/5Pseq/deduplicated/pt221-lo_Aligned.sortedByCoord.dedup.bam \
    ~/data/tss/human/5Pseq/multiqc/multiqc_report.html \
    ~/data/tss/human/5Pseq/alignment_qc/multiqc_report.html |  dot -Tpng > dag/alignment_dag.png
