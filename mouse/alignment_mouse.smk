# how to run from ./:
# snakemake -s alignment.smk --jobs 5000 --latency-wait 120 --cluster-config config/cluster.json --cluster 'qsub -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e  {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_alignment.yml"

SAMPLE=['mESC1', 'mESC2', 'mESC3']
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss/mouse'

rule all:
    input:
        expand("{pdir}/1_raw_data/5Pseq_processed/fastq/{sample}_{replicate}-htseq-qc.pdf",
            replicate=[1,2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{genomedir}/STARINDEX/Genome",
            genomedir=config['genome']),
        expand("{pdir}/2_alignments/{sample}_{replicate}_Aligned.out.sam",
            replicate=[1,2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/2_alignments/STAR_summary.csv",
            pdir=DIRECTORY)

rule qc:
    input:
        "{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_{replicate}.fastq",
    output:
        "{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_{replicate}-htseq-qc.pdf",
    shell:
        "htseq-qa --type=fastq --outfile={output} {input}"

rule generate_genome:
    input:
        genome="{genomedir}/NCBIM37.genome.fa",
        gtf="{genomedir}/gencode.vM1.annotation.gtf"
    output:
        genome="{genomedir}/STARINDEX/Genome"
    params:
        length=49,
        thread=8,
    shell:
        """
        STAR --runThreadN {params.thread} \
        --runMode genomeGenerate \
        --genomeDir {wildcards.genomedir}/STARINDEX \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {wildcards.genomedir}/gencode.vM1.annotation.gtf \
        --sjdbOverhang {params.length}
        """

rule align:
    input:
        read="{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_{replicate}.fastq",
    output:
        sam="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.sam",
        out="{dir}/2_alignments/{sample}_{replicate}_Log.final.out",
    params:
        thread=4,
        genomedir=config['genome']
    shell:
        """
        STAR --runThreadN {params.thread} \
        --runMode alignReads \
        --genomeDir {params.genomedir}/STARINDEX \
        --readFilesIn {input.read} \
        --outReadsUnmapped Fastq \
        --outFileNamePrefix {wildcards.dir}/2_alignments/{wildcards.sample}_{wildcards.replicate} \
        """

rule overview_STAR_results:
    input:
        positions=expand("{{dir}}/2_alignments/{sample}_{replicate}_Log.final.out",
            sample=SAMPLE,
            replicate=[1,2]),
    output:
        csv="{dir}/2_alignments/STAR_summary.csv",
        pdf="{dir}/2_alignments/STAR_summary.pdf",
    shell:
        """
        Rscript alignment/STAR-alignment-results.r \
            --directory {wildcards.dir}/2_alignments \
            --suffix _Log.final.out \
            --verbose
        """
