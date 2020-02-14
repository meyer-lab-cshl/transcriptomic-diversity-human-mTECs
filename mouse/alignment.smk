# how to run from ./:
# snakemake -s alignment.smk --jobs 5000 --latency-wait 120 --cluster-config config/cluster.json --cluster 'qsub -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e  {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_alignment.yml"

SAMPLES=['mESC1', 'mESC2', 'mESC3']
REPLICATES=[1,2]
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss/mouse/5Pseq'

rule all:
    input:
        expand("{pdir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam.bai",
            replicate=REPLICATES,
            pdir=DIRECTORY,
            sample=SAMPLES),
        expand("{pdir}/multiqc/multiqc_report.html",
            pdir=DIRECTORY),
        expand("{pdir}/alignment_qc/multiqc_report.html",
            pdir=DIRECTORY)


rule fastqc_raw:
    input:
        samples="{dir}/fastq/{sample}_{replicate}.fastq"
    output:
        "{dir}/fastqc_raw/{sample}_{replicate}_fastqc.html"
    threads:
        2
    shell:
        """
        module load FastQC/0.11.8-Java-1.8
        fastqc \
            -t {threads} \
            -o "{wildcards.dir}/fastqc_raw" \
            {input.samples:q}
        """

rule fastq_screen:
    input:
        reads="{dir}/fastq/{sample}_{replicate}.fastq",
    output:
        screen="{dir}/fastq_screen/{sample}_{replicate}_screen.txt",
        filtered="{dir}/fastq_filtered/{sample}_{replicate}.fastq",
    threads:
        2
    shell:
        """
        # Filter: mouse reads, second position in filter string
        # unique+multi mappers (3)
        fastq_screen \
            --force \
            --threads {threads} \
            --tag \
            --filter 03000000000000 \
            --outdir {wildcards.dir}/fastq_screen \
            {input.reads:q}
        mv {wildcards.dir}/fastq_screen/{wildcards.sample}_{wildcards.replicate}.tagged_filter.fastq \
           {output.filtered}
        """

rule fastqc_filtered:
    input:
        samples="{dir}/fastq_filtered/{sample}_{replicate}.fastq"
    output:
        "{dir}/fastqc_filtered/{sample}_{replicate}_fastqc.html"
    threads:
        2
    shell:
        """
        module load FastQC/0.11.8-Java-1.8
        fastqc \
            -t {threads} \
            -o "{wildcards.dir}/fastqc_filtered" \
            {input.samples:q}
        """

rule multiqc_raw:
    input:
        expand("{{dir}}/fastq_screen/{sample}_{replicate}_screen.txt",
            sample=SAMPLES,
            replicate=REPLICATES),
        expand("{{dir}}/fastqc_filtered/{sample}_{replicate}_fastqc.html",
            sample=SAMPLES,
            replicate=REPLICATES)
    output:
        "{dir}/multiqc/multiqc_report.html"
    shell:
        """
        multiqc \
            --force \
            --export \
            --outdir {wildcards.dir}/multiqc \
            {wildcards.dir}/fastq_screen/ \
            {wildcards.dir}/fastqc_filtered/
        """

rule align:
    input:
        read="{dir}/fastq_filtered/{sample}_{replicate}.fastq",
        genome=expand("{genomedir}/GRCm38/STARINDEX/Genome",
            genomedir=config['genome'])
    output:
        bam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam",
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
        --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix {wildcards.dir}/alignments/{wildcards.sample}_{wildcards.replicate}_ \
        """

rule index:
    input:
        bam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam",
    output:
        index="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam.bai",
    shell:
        """
        samtools index {input.bam}
        """

rule rnaseq_metrics_raw:
    input:
        bam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam",
        refflat=config['refflat']
    output:
        metrics="{dir}/metrics/{sample}_{replicate}.rna_metrics"
    shell:
        """
        module load picard/2.18.20
        java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
            VALIDATION_STRINGENCY=LENIENT \
            REF_FLAT={input.refflat:q} \
            STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
            ASSUME_SORTED=false \
            INPUT={input.bam:q} \
            OUTPUT={output.metrics:q}
        """

rule multiqc_aligned:
    input:
        expand("{{dir}}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam",
            replicate=REPLICATES,
            sample=SAMPLES),
        expand("{{dir}}/metrics/{sample}_{replicate}.rna_metrics",
            replicate=REPLICATES,
            sample=SAMPLES),
    output:
        report="{dir}/alignment_qc/multiqc_report.html"
    shell:
        """
        multiqc \
            --force \
            --export \
            --outdir {wildcards.dir:q}/alignment_qc \
               {wildcards.dir}/alignments {wildcards.dir}/metrics
        """

