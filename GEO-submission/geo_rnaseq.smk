
SAMPLES=['pt214_hi', 'pt214_lo', 'pt221_hi', 'pt221_lo', 'pt226_lo', 'pt226_hi']


rule all:
    input:
        expand("{pdir}/GEO/{sample}_Aligned.sortedByCoord.out.bam",
            pdir="/grid/meyer/home/jacarter/TSS_RNAseq",
            sample=SAMPLES),


rule bamboozle:
    input:
        bam="{pdir}/Aligned/{sample}_Aligned.sortedByCoord.out.bam",
        fasta="/grid/meyer/home/common/public/annotations/genome/human/GRCh38/human.GRCh38.genome.fa",
        faidx="/grid/meyer/home/common/public/annotations/genome/human/GRCh38/human.GRCh38.genome.fa.fai"
    output:
        bam="{pdir}/GEO/{sample}_Aligned.sortedByCoord.out.bam",
    threads: 4
    shell:
        """
        BAMboozle \
            --bam {input.bam} \
            --out {output} \
            --fa {input.fasta} \
            --p {threads} \
            --strict \
            --keepsecondary \
            --keepunmapped
        """

rule convertfastq:
    input:
        bam="{pdir}/GEO/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        r1="{pdir}/GEO/{sample}_r1_fastq.gz",
        r2="{pdir}/GEO/{sample}_r2_fastq.gz",
        ambiguous="{pdir}/GEO/{sample}_ambiguous_fastq.gz",
        singleton="{pdir}/GEO/{sample}_singleton_fastq.gz",
    shell:
        """
        # sort paired read alignment .bam file (sort by name -n)
        samtools sort -n {input.bam} -o {input.bam}.tmp

        # save fastq reads in separate R1 and R2 files
        samtools fastq -@ 8 {input.bam}.tmp \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 {output.ambiguous} \
            -s {output.singleton} \
            -n

