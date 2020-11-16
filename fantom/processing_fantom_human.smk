# how to run from ./:
# snakemake -p -s processing_fantom_human.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json
# --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

#configfile: "config/config_processing_alignment.yml"

DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss/human/fantom'
ANNOTATION='/sonas-hs/meyer/hpc/home/hmeyer/data/common/public/annotations/genome/liftover'
SAMPLES, = glob_wildcards("/sonas-hs/meyer/hpc/home/hmeyer/data/tss/human/fantom/bed/GRCh37/{sample}.hg19.ctss.bed.gz")

rule all:
    input:
        expand("{pdir}/bed/GRCh38/{sample}.bed",
            sample=SAMPLES,
            pdir=DIRECTORY),
        expand("{pdir}/tss/combined/all_tissues.positions.csv",
            pdir=DIRECTORY),
        expand("{pdir}/tss/summary/all_tissues_consensus_sites.csv",
            pdir=DIRECTORY),
        expand("{pdir}/tss/isoforms/{sample}.isoforms.csv",
            sample=SAMPLES,
            pdir=DIRECTORY)
        #"{dir}/combined/fantom_liftover_human_mouse.pdf",


rule liftover_counts:
    input:
        positions="{dir}/bed/GRCh37/{sample}.hg19.ctss.bed.gz",
        chain=expand("{annodir}/hg19ToHg38.over.chain",
            annodir=ANNOTATION),
    output:
        liftover="{dir}/bed/GRCh38/{sample}.bed",
        unmapped="{dir}/bed/GRCh38/{sample}.unmapped.bed",
        count_unmapped="{dir}/bed/GRCh38/{sample}.count.unmapped.txt",
    shell:
        """
        liftOver {input.positions} {input.chain} \
            {output.liftover}.tmp {output.unmapped}.tmp
        awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2".."$3","$6,$5,$6}}' {output.liftover}.tmp \
            > {output.liftover}
        awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2".."$3","$6,$5,$6}}' {output.unmapped}.tmp \
            > {output.unmapped}
        echo "{wildcards.sample}" > {output.count_unmapped}.tmp
        wc -l {output.liftover}  | cut -d " " -f 1 \
            > {output.count_unmapped}.mapped.tmp
        grep "#" -v {output.unmapped} | wc -l | \
            paste {output.count_unmapped}.tmp {output.count_unmapped}.mapped.tmp - \
            > {output.count_unmapped}
        rm {output.liftover}.tmp {output.unmapped}.tmp \
            {output.count_unmapped}.tmp {output.count_unmapped}.mapped.tmp
        """

rule summarise_liftover:
    input:
        expand("{{dir}}/GRCh38/{sample}.count.unmapped.txt",
            sample=SAMPLES)
    output:
        "{dir}/bed/GRCh38/all_tissues.count.unmapped.txt",
    shell:
        """
        cat {input} > {output}
        """

rule process_counts:
    input:
        counts="{dir}/bed/GRCh38/{sample}.bed",
    output:
        bedgraph="{dir}/tss/bedgraphs/{sample}.bedgraph",
        counts="{dir}/tss/summary/{sample}.summary.counts.csv",
        positions="{dir}/tss/raw_positions/{sample}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --species human \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/tss \
            --sample {wildcards.sample}
        """

rule combine_counts:
    input:
        positions=expand("{{dir}}/tss/raw_positions/{sample}.positions.csv",
            sample=SAMPLES)
    output:
        positions="{dir}/tss/combined/all_tissues.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/tss/raw_positions \
            --odir {wildcards.dir}/tss/combined \
            --suffix .positions.csv \
            --fantom \
            --verbose
        """


rule visualise_liftover:
    input:
        human="{dir}/human/fantom/bed/GRCh38/all_tissues.count.unmapped.txt",
        mouse="{dir}/mouse/fantom/bed/GRCm38/all_mESC_C46.count.unmapped.txt",
    output:
        "{dir}/combined/fantom_liftover_human_mouse.pdf"
    shell:
        """
        Rscript liftover_summary.R \
            --dir {wildcards.dir}/combined \
            --human {input.human} \
            --mouse {input.mouse}
        """


rule consensus_sites:
    input:
        summary=expand("{{dir}}/tss/summary/{sample}.summary.counts.csv",
            sample=SAMPLES),
    output:
        collapsed="{dir}/tss/summary/all_tissues_consensus_sites.csv",
        summary="{dir}/tss/summary/all_tissues_consensus_sites_summary.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_consensus.pl\
            --dir {wildcards.dir}/tss/summary \
            --collapsed {output.collapsed} \
            --summarised {output.summary} \
            --maxdist 12\
            --minbc 10\
            --suffix .summary.counts.csv \
        """

rule map_sample_sites:
    input:
        summary="{dir}/tss/summary/{sample}.summary.counts.csv",
        collapsed="{dir}/tss/summary/all_tissues_consensus_sites.csv",
    output:
        isoform="{dir}/tss/isoforms/{sample}.isoforms.csv",
        gene="{dir}/tss/genes/{sample}.genes.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_sample.pl\
            --insummary {input.summary} \
            --iniso {input.collapsed} \
            --outiso {output.isoform} \
            --outgene {output.gene} \
        """
