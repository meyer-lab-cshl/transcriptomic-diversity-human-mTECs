snakemake -s alignment_mouse.smk \
    --dag ~/data/tss/mouse/2_alignments/mESC1_1_Aligned.out.sam \
    ~/data/tss/mouse/2_alignments/STAR_summary.pdf \
    ~/data/common/public/annotations/mouse/genome/NCBI37_mm9/STARINDEX/Genome|  dot -Tpdf > alignment_mouse_dag.pdf
