
SAMPLES=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
         'pt221-hi', 'pt221-lo', 'pt226-lo', 'pt226-hi']


rule all:
    input:
        expand("{pdir}/GEO/{sample}_Aligned.sortedByCoord.out.bam",
            pdir="/grid/meyer/home/hmeyer/data/tss/human/5Pseq",
            sample=SAMPLES),


rule bamboozle:
    input:
        bam="{pdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
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


