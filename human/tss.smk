# how to run from ./:
# snakemake -s tss.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
        'pt221-hi', 'pt221-lo', 'pt226-hi', 'pt226-lo']
DIRECTORY='/mnt/grid/meyer/hpc/home/data/hmeyer/data/tss/human/5Pseq'

rule all:
    input:
        expand("{pdir}/tss/summary/all_samples_consensus_sites.csv",
            pdir=DIRECTORY),
        expand("{pdir}/tss/isoforms/{sample}.isoforms.csv",
            sample=SAMPLE,
            pdir=DIRECTORY)

rule consensus_sites:
    input:
        summary=expand("{{dir}}/tss/summary/{sample}_fwd.summary.counts.csv",
            sample=SAMPLE),
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
            --suffix _fwd.summary.counts.csv \
        """

rule map_sample_sites:
    input:
        summary="{dir}/tss/summary/{sample}_fwd.summary.counts.csv",
        collapsed="{dir}/tss/summary/all_samples_consensus_sites.csv",
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
