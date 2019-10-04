# how to run from ./:
# snakemake -s processing_alignment.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
        'pt221-hi', 'pt221-lo', 'pt226-hi', 'pt226-lo']
DIRECTORY='/mnt/grid/meyer/hpc/home/data/tss'

rule all:
    input:
        expand("{pdir}/combined/3_tss_data/summary/all_mTECs_consensus_sites.csv",
            pdir=DIRECTORY)

rule consensus_sites:
    input:
        summary=expand("{{dir}}/human/3_tss_data/summary/{sample}_fwd.summary.counts",
            sample=SAMPLE),
    output:
        collapsed="{dir}/combined/3_tss_data/summary/all_mTECs_consensus_sites.csv",
        summary="{dir}/combined/3_tss_data/summary/all_mTECs_consensus_sites_summary.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_consensus.pl\
            --dir {wildcards.dir}/human/3_tss_data/raw_positions \
            --collapsed {output.collapsed} \
            --summarised {output.summary} \
            --maxdist 12\
            --minbc 10\
            --suffix _fwd.summary.counts \
        """

rule map_sample_sites:
    input:
        summary="{dir}/human/3_tss_data/summary/{sample}_fwd.summary.counts",
        collapsed="{dir}/combined/3_tss_data/summary/all_mTECs_consensus_sites.csv",
    output:
        isoform="{dir}/human/3_tss_data/isoforms/{sample}.isoforms.csv",
        gene="{dir}/human/3_tss_data/genes/{sample}.genes.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_sample.pl\
            --iniso {input.collapsed} \
            --outiso {output.isoform} \
            --outgene {output.gene} \
        """
