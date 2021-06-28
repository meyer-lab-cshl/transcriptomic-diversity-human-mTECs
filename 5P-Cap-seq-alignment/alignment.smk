# how to run from ./:
# snakemake -np -s alignment.smk --jobs 5000 --latency-wait 180 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_alignment.yml"

DIRECTORY="/sonas-hs/meyer/hpc/home/hmeyer/data/tss/human/5Pseq"
#SAMPLES=['pt87-hi']
SAMPLES=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
         'pt221-hi', 'pt221-lo', 'pt226-lo', 'pt226-hi']

rule all:
    input:
        expand("{pdir}/multiqc/multiqc_report.html",
            pdir=DIRECTORY),
        expand("{pdir}/deduplicated/{sample}_Aligned.sortedByCoord.{suffix}.bam",
            suffix=['dedup'],
            pdir=DIRECTORY,
            sample=SAMPLES),
        expand("{pdir}/alignment_qc/multiqc_report.html",
            pdir=DIRECTORY),

rule demultiplex:
    input:
        read1="{readsdir}/raw_data/mTEC-multiplexed_1.txt.gz",
        read2="{readsdir}/raw_data/mTEC-multiplexed_2.txt.gz",
        barcodes="{readsdir}/raw_data/mTEC_5Seq_2014-03-18.design.fasta",
    output:
        pt87_lo_1="{readsdir}/fastq/pt87-lo_1.fastq",
        pt212_lo_1="{readsdir}/fastq/pt212-lo_1.fastq",
        pt214_lo_1="{readsdir}/fastq/pt214-lo_1.fastq",
        pt87_hi_1="{readsdir}/fastq/pt87-hi_1.fastq",
        pt212_hi_1="{readsdir}/fastq/pt212-hi_1.fastq",
        pt214_hi_1="{readsdir}/fastq/pt214-hi_1.fastq",
        pt87_lo_2="{readsdir}/fastq/pt87-lo_2.fastq",
        pt212_lo_2="{readsdir}/fastq/pt212-lo_2.fastq",
        pt214_lo_2="{readsdir}/fastq/pt214-lo_2.fastq",
        pt87_hi_2="{readsdir}/fastq/pt87-hi_2.fastq",
        pt212_hi_2="{readsdir}/fastq/pt212-hi_2.fastq",
        pt214_hi_2="{readsdir}/fastq/pt214-hi_2.fastq",
    params:
    shell:
        """
        #input order swapped as barcode on reverse read
        cutadapt \
            -e 0.2 \
            --no-indels \
            -g file:{input.barcodes} \
            -o {wildcards.readsdir}/fastq/{{name}}_2.fastq \
            -p {wildcards.readsdir}/fastq/{{name}}_1.fastq \
            {input.read2} {input.read1}
        """

rule format:
    input:
        read1="{readsdir}/raw_data/{sample}_1.txt.gz",
        read2="{readsdir}/raw_data/{sample}_2.txt.gz"
    output:
        read1="{readsdir}/fastq/{sample}_1.fastq",
        read2="{readsdir}/fastq/{sample}_2.fastq",
    wildcard_constraints:
        sample="pt221-hi|pt221-lo|pt226-hi|pt226-lo"
    shell:
        """
        zcat {input.read1} > {output.read1}
        zcat {input.read2} > {output.read2}
        """

rule fastqc_raw:
    input:
        samples="{readsdir}/fastq/{sample}_{read}.fastq"
    output:
        "{readsdir}/fastqc_raw/{sample}_{read}_fastqc.html"
    threads:
        2
    shell:
        """
        module load FastQC/0.11.8-Java-1.8
        fastqc \
            -t {threads} \
            -o "{wildcards.readsdir}/fastqc_raw" \
            {input.samples:q}
        """

rule umi:
    input:
        read1="{readsdir}/fastq/{sample}_1.fastq",
        read2="{readsdir}/fastq/{sample}_2.fastq",
    output:
        read1="{readsdir}/fastq_trimmed/{sample}_1.fastq",
        read2="{readsdir}/fastq_trimmed/{sample}_2.fastq",
    params:
        extract="NNNNNNNN",
    shell:
        """
        umi_tools extract \
            --extract-method=string \
            --bc-pattern={params.extract} \
            -I {input.read1} \
            --read2-in={input.read2} \
            --stdout={output.read1} \
            --read2-out={output.read2} \
            -L {wildcards.readsdir}/extract.log
        """

rule fastqc_trimmed:
    input:
        samples="{readsdir}/fastq_trimmed/{sample}_{read}.fastq"
    output:
        "{readsdir}/fastqc_trimmed/{sample}_{read}_fastqc.html"
    threads:
        2
    shell:
        """
        module load FastQC/0.11.8-Java-1.8
        fastqc \
            -t {threads} \
            -o "{wildcards.readsdir}/fastqc_trimmed" \
            {input.samples:q}
        """

rule fastq_screen:
    input:
        reads="{readsdir}/fastq_trimmed/{sample}_{read}.fastq",
    output:
        screen="{readsdir}/fastq_screen/{sample}_{read}.tagged_filter.fastq",
        filtered="{readsdir}/fastq_filtered/{sample}_{read}.fastq",
    threads:
        2
    shell:
        """
        # Filter: human reads, first position in filter string
        # unique+multi mappers (3)
        fastq_screen \
            --force \
            --threads {threads} \
            --tag \
            --filter 30000000000000 \
            --outdir {wildcards.readsdir}/fastq_screen \
            {input.reads:q}
        mv {output.screen} {output.filtered}
        """

rule fastq_pair:
    input:
        read1="{readsdir}/fastq_filtered/{sample}_1.fastq",
        read2="{readsdir}/fastq_filtered/{sample}_2.fastq",
    output:
        read1="{readsdir}/fastq_filtered/{sample}_1.paired.fastq",
        read2="{readsdir}/fastq_filtered/{sample}_2.paired.fastq",
    params:
        t=5000000
    threads:
        2
    shell:
        """
        fastq_pair \
            -t {params.t} \
            -p {input.read1} {input.read2}
        mv {input.read1}.paired.fq {output.read1}
        mv {input.read2}.paired.fq {output.read2}
        """

rule fastqc_filtered:
    input:
        samples="{readsdir}/fastq_filtered/{sample}_{read}.paired.fastq"
    output:
        "{readsdir}/fastqc_filtered/{sample}_{read}.paired_fastqc.html"
    threads:
        2
    shell:
        """
        module load FastQC/0.11.8-Java-1.8
        fastqc \
            -t {threads} \
            -o "{wildcards.readsdir}/fastqc_filtered" \
            {input.samples:q}
        """

rule multiqc_raw:
    input:
        expand("{{readsdir}}/fastqc_trimmed/{sample}_{read}_fastqc.html",
            sample=SAMPLES,
            read=[1,2]),
        expand("{{readsdir}}/fastq_screen/{sample}_{read}_screen.txt",
            sample=SAMPLES,
            read=[1,2]),
        expand("{{readsdir}}/fastqc_filtered/{sample}_{read}.paired_fastqc.html",
            sample=SAMPLES,
            read=[1,2])
    output:
        "{readsdir}/multiqc/multiqc_report.html"
    shell:
        """
        multiqc \
            --force \
            --export \
            --outdir {wildcards.readsdir}/multiqc \
            {wildcards.readsdir}/fastqc_trimmed/  \
            {wildcards.readsdir}/fastq_screen/ \
            {wildcards.readsdir}/fastqc_filtered/
        """

rule align:
    input:
        read1="{readsdir}/fastq_filtered/{sample}_1.paired.fastq",
        read2="{readsdir}/fastq_filtered/{sample}_2.paired.fastq",
        genome=expand("{genomedir}/STARINDEX/Genome",
            genomedir=config['genomehuman'])
    output:
        "{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        "{readsdir}/alignments/{sample}_Log.final.out",
    params:
        thread=4,
        genomedir=config['genomehuman']
    shell:
        """
        STAR --runThreadN {params.thread} \
            --runMode alignReads \
            --genomeDir {params.genomedir}/STARINDEX \
            --readFilesIn {input.read1} {input.read2} \
            --outReadsUnmapped Fastq \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outFileNamePrefix {wildcards.readsdir}/alignments/{wildcards.sample}_ \
        """

rule index:
    input:
        bam="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        index="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai",
    shell:
        """
        samtools index {input.bam}
        """

rule umi_group:
    input:
        bam="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        index="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        group="{readsdir}/grouped/{sample}_Aligned.sortedByCoord.group.tsv",
        bam="{readsdir}/grouped/{sample}_Aligned.sortedByCoord.group.bam",
    shell:
        """
        umi_tools group \
            -I {input.bam} \
            --paired \
            --group-out={output.group} \
            --output-bam \
            -S {output.bam}
        """

rule dedup:
    input:
        bam="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        index="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        bam="{readsdir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.bam",
        #stats="{readsdir}/deduplicated/{sample}_Aligned.sortedByCoord.stats_per_umi.tsv"
    shell:
        """
        umi_tools dedup \
            -I {input.bam} \
            -L {wildcards.readsdir}/deduplicated/extract.log \
            -S {output.bam} \
            --paired \
            --output-stats={wildcards.readsdir}/deduplicated/{wildcards.sample}_Aligned.sortedByCoord.stats_
        """



rule rnaseq_metrics_raw:
    input:
        bam="{readsdir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        refflat=config['refflat']
    output:
        metrics="{readsdir}/metrics/{sample}.rna_metrics"
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

rule rnaseq_metrics_dedup:
    input:
        bam="{readsdir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.bam",
        refflat=config['refflat']
    output:
        metrics="{readsdir}/metrics_deduplicated/{sample}.deduplicated.rna_metrics"
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
        expand("{{readsdir}}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
            sample=SAMPLES),
        expand("{{readsdir}}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.bam",
            sample=SAMPLES),
        expand("{{readsdir}}/metrics/{sample}.rna_metrics",
            sample=SAMPLES),
        expand("{{readsdir}}/metrics_deduplicated/{sample}.deduplicated.rna_metrics",
            sample=SAMPLES),
    output:
        report="{readsdir}/alignment_qc/multiqc_report.html"
    shell:
        """
        multiqc \
            --force \
            --export \
            --outdir {wildcards.readsdir:q}/alignment_qc \
               {wildcards.readsdir}/deduplicated {wildcards.readsdir}/metrics_deduplicated
        """

#{wildcards.readsdir}/alignments {wildcards.readsdir}/metrics \
