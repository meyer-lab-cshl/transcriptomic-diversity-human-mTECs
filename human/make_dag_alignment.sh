snakemake -s alignment_human.smk \
    --dag ~/data/tss/human/1_raw_data/5Pseq_processed/fastq/pt221-lo_1.fastq \
    ~/data/tss/human/2_alignments/pt221-lo_Aligned.out.sam \
    ~/data/tss/human/2_alignments/STAR_summary.pdf \
    ~/data/common/public/annotations/human/genome/STARINDEX/Genome|  dot -Tpdf > alignment_human_dag.pdf
