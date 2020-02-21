# how to run from ./:
# snakemake -s tss.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

REPLICATES=[1,2]
SAMPLES=['mESC1', 'mESC2', 'mESC3']
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss/mouse/5Pseq'

rule all:
    input:
        expand("{pdir}/tss/summary/all_samples_consensus_sites.csv",
            pdir=DIRECTORY),
        expand("{pdir}/tss/isoforms/{sample}_{replicate}.isoforms.csv",
            sample=SAMPLES,
            replicate=REPLICATES,
            pdir=DIRECTORY)

rule consensus_sites:
    input:
        summary=expand("{{dir}}/tss/summary/{sample}_{replicate}.summary.counts.csv",
            replicate=REPLICATES,
            sample=SAMPLES),
    output:
        collapsed="{dir}/tss/summary/all_samples_consensus_sites.csv",
        summary="{dir}/tss/summary/all_samples_consensus_sites_summary.csv",
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
        summary="{dir}/tss/summary/{sample}_{replicate}.summary.counts.csv",
        collapsed="{dir}/tss/summary/all_samples_consensus_sites.csv",
    output:
        isoform="{dir}/tss/isoforms/{sample}_{replicate}.isoforms.csv",
        gene="{dir}/tss/genes/{sample}_{replicate}.genes.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_sample.pl\
            --insummary {input.summary} \
            --iniso {input.collapsed} \
            --outiso {output.isoform} \
            --outgene {output.gene} \
        """
