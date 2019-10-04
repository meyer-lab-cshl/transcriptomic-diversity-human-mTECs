# how to run from ./:
# snakemake -s processing_alignment.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['mESC1', 'mESC2', 'mESC3']
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss'

rule all:
    input:
        expand("{pdir}/combined/3_tss_data/raw_positions/single_samples/{sample}.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/combined/3_tss_data/raw_positions/all_tissues.hg19.fantom.positions.csv",
            pdir=DIRECTORY)

rule process_counts_mouse:
    input:
        counts="{dir}/fantom/mouse/{sample}.mm9.fantom.bed",
    output:
        bedgraph="{dir}/combined/3_tss_data/bedgraphs/single_samples/{sample}.mm9.fantom.bedgraph",
        counts="{dir}/combined/3_tss_data/summary/single_samples/{sample}.mm9.fantom.summary.counts.csv",
        positions="{dir}/combined/3_tss_data/raw_positions/single_samples/{sample}.mm9.fantom.positions.csv",
    shell:
        """
        Rscript processing/process_counts.r \
            --species mouse \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/combined/3_tss_data \
            --sample {wildcards.sample}
        """

rule combine_counts:
    input:
        positions=expand("{{dir}}/3_tss_data/raw_positions/single_samples/{sample}.mm9.fantom.positions.csv",
            sample=SAMPLE,
            replicate=[1,2]),
    output:
        positions="{dir}/3_tss_data/raw_positions/all_tissues.mm9.fantom.positions.csv",
    shell:
        """
        Rscript processing_alignment/combine_counts.r \
            --directory {wildcards.dir}/3_tss_data/raw_positions \
            --suffix .mm9.fantom.positions.csv \
            --ofile {output.positions} \
            --verbose
        """

rule process_counts_human:
    input:
        counts="{dir}/fantom/mouse/{sample}.hg19.fantom.bed",
    output:
        bedgraph="{dir}/combined/3_tss_data/bedgraphs/single_samples/{sample}.hg19.fantom.bedgraph",
        counts="{dir}/combined/3_tss_data/summary/single_samples/{sample}.hg19.fantom.summary.counts.csv",
        positions="{dir}/combined/3_tss_data/raw_positions/single_samples/{sample}.hg19.fantom.positions.csv",
    shell:
        """
        Rscript processing/process_counts.r \
            --species mouse \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/combined/3_tss_data \
            --sample {wildcards.sample}
        """

rule combine_counts_human:
    input:
        positions=expand("{{dir}}/3_tss_data/raw_positions/single_samples/{sample}.hg19.fantom.positions.csv",
            sample=SAMPLE,
            replicate=[1,2]),
    output:
        positions="{dir}/3_tss_data/raw_positions/all_tissues.hg19.fantom.positions.csv",
    shell:
        """
        Rscript processing_alignment/combine_counts.r \
            --directory {wildcards.dir}/3_tss_data/raw_positions \
            --suffix .hg19.fantom.positions.csv \
            --ofile {output.positions} \
            --verbose
        """
