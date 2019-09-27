# how to run from ./:
# snakemake -s alignment_human_thymus.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_alignment.yml"

SAMPLE=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
        'pt221-hi', 'pt221-lo', 'pt226-hi', 'pt226-lo']
DIRECTORY='/mnt/grid/meyer/hpc/home/data/common/tss'

rule all:
    input:
        expand("{pdir}/1_raw_data/5Pseq_processed/fastq/{sample}_{read}-htseq-qc.pdf",
            pdir=DIRECTORY,
            read=[1,2],
            sample=SAMPLE),
        expand("{genomedir}/STARINDEX/Genome",
            genomedir=config['genomehuman']),
        expand("{pdir}/2_alignments/{sample}_{suffix}",
            suffix=["Aligned.out.sam", "Log.final.out"],
            pdir=DIRECTORY,
            sample=SAMPLE),

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.demultiplexing.get(**wildcards).output[0]
    return expand("{dir}/1_raw_data/5Pseq_processed/{file}/{samples}_1.fastq",
                  samples=glob_wildcards(os.path.join(checkpoint_output,
                                                  "{samples}_1.fastq")).samples)


rule demultiplexing_and_trim:
    input:
        read1="{dir}/1_raw_data/5Pseq/mTEC-multiplexed_1.txt.gz",
        read2="{dir}/1_raw_data/5Pseq/mTEC-multiplexed_2.txt.gz",
        design="{dir}/1_raw_data/5Pseq/mTEC-multiplexed.design",
    output:
        pt87_lo_1="{dir}/1_raw_data/5Pseq_processed/fastq/pt87-lo_1.fastq",
        pt212_lo_1="{dir}/1_raw_data/5Pseq_processed/fastq/pt212-lo_1.fastq",
        pt214_lo_1="{dir}/1_raw_data/5Pseq_processed/fastq/pt214-lo_1.fastq",
        pt87_hi_1="{dir}/1_raw_data/5Pseq_processed/fastq/pt87-hi_1.fastq",
        pt212_hi_1="{dir}/1_raw_data/5Pseq_processed/fastq/pt212-hi_1.fastq",
        pt214_hi_1="{dir}/1_raw_data/5Pseq_processed/fastq/pt214-hi_1.fastq",
        pt87_lo_2="{dir}/1_raw_data/5Pseq_processed/fastq/pt87-lo_2.fastq",
        pt212_lo_2="{dir}/1_raw_data/5Pseq_processed/fastq/pt212-lo_2.fastq",
        pt214_lo_2="{dir}/1_raw_data/5Pseq_processed/fastq/pt214-lo_2.fastq",
        pt87_hi_2="{dir}/1_raw_data/5Pseq_processed/fastq/pt87-hi_2.fastq",
        pt212_hi_2="{dir}/1_raw_data/5Pseq_processed/fastq/pt212-hi_2.fastq",
        pt214_hi_2="{dir}/1_raw_data/5Pseq_processed/fastq/pt214-hi_2.fastq",
    params:
        sample_bc=6,
        molecular_bc=8,
        info=lambda wildcards: "{}/1_raw_data/5Pseq_processed/fastq/mTEC-multiplexed.info".format(wildcards.dir)
    shell:
        """
        python alignment/sort-and-trim-mbc.py \
            --forward {input.read1} \
            --reverse {input.read2} \
            --design {input.design} \
            --sample-barcode {params.sample_bc} \
            --molecular-barcode {params.molecular_bc} \
            --trim-molecular-barcode \
            --multiplexed \
            --info {params.info}
        """

rule trim:
    input:
        read1="{dir}/1_raw_data/5Pseq/{sample}_1.txt.gz",
        read2="{dir}/1_raw_data/5Pseq/{sample}_2.txt.gz",
        design="{dir}/1_raw_data/5Pseq/{sample}.design",
    output:
        pt_1="{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_1.fastq",
        pt_2="{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_2.fastq",
        info="{dir}/1_raw_data/5Pseq_processed/fastq/{sample}.info",
    wildcard_constraints:
        sample="pt221-hi|pt221-lo|pt226-hi|pt226-lo"
    params:
        sample_bc=6,
        molecular_bc=8,
    shell:
        """
        python alignment/sort-and-trim-mbc.py \
            --forward {input.read1} \
            --reverse {input.read2} \
            --design {input.design} \
            --sample-barcode {params.sample_bc} \
            --molecular-barcode {params.molecular_bc} \
            --trim-molecular-barcode \
            --info {output.info}
        """
rule qc:
    input:
        "{dir}/1_raw_data/5Pseq_processed/fastq/{sample}-{read}.fastq",
    output:
        "{dir}/1_raw_data/5Pseq_processed/fastq/{sample}-{read}-htseq-qc.pdf",
    shell:
        "htseq-qa --type=fastq --outfile={output} {input}"

rule add_spikeins_genome:
    input:
        genome="{genomedir}/GRCh38.primary_assembly.genome.fa",
        spikeins=expand("{{genomedir}}/pGIBS-{type}.fasta",
            type=config['spikeins'])
    output:
        genome="{genomedir}/GRCh38.primary_assembly.genome.spikeins.fa",
    shell:
        "cat {input.genome} {input.spikeins} > {output.genome}"

rule generate_genome:
    input:
        genome="{genomedir}/GRCh38.primary_assembly.genome.spikeins.fa",
        gtf="{genomedir}/gencode.v31.primary_assembly.annotation.gtf"
    output:
        genome="{genomedir}/STARINDEX/Genome"
    params:
        length=49,
        thread=4,
    shell:
        """
        STAR --runThreadN {params.thread} \
        --runMode genomeGenerate \
        --genomeDir {wildcards.genomedir}/STARINDEX \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.length}
        """

rule align:
    input:
        read1="{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_1.fastq",
        read2="{dir}/1_raw_data/5Pseq_processed/fastq/{sample}_2.fastq",
    output:
        "{dir}/2_alignments/{sample}_Aligned.out.sam",
        "{dir}/2_alignments/{sample}_Log.final.out",
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
        --outFileNamePrefix {wildcards.dir}/2_alignments/{wildcards.sample}_ \
        """
