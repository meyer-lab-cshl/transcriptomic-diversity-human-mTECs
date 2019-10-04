snakemake -s alignment_human.smk \
    --dag /mnt/grid/meyer/hpc/home/data/common/tss/1_raw_data/5Pseq_processed/fastq/pt221-lo_1.fastq \
    /mnt/grid/meyer/hpc/home/data/common/tss/2_alignments/pt221-lo_Aligned.out.sam \
    /mnt/grid/meyer/hpc/home/data/common/tss/2_alignments/STAR_summary.pdf \
    /mnt/grid/meyer/hpc/home/data/common/public/annotations/human/genome/STARINDEX/Genome|  dot -Tpdf > alignment_human_dag.pdf
