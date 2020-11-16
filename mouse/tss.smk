# how to run from ./:
# snakemake -s tss.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

REPLICATES=[1,2]
SAMPLES=['mESC1', 'mESC2', 'mESC3']
FANTOM=['ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1',
'ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep2.CNhs14109.14362-155I6']
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss'

rule all:
    input:
        expand("{pdir}/combined/mouse_combined_consensus_sites.csv",
            pdir=DIRECTORY),
        expand("{pdir}/mouse/5Pseq/tss/isoforms/{sample}_{replicate}.isoforms.csv",
            sample=SAMPLES,
            replicate=REPLICATES,
            pdir=DIRECTORY),
        expand("{pdir}/mouse/fantom/tss/isoforms/{fantom}.isoforms.csv",
            fantom=FANTOM,
            pdir=DIRECTORY),
        expand("{pdir}/combined/tss/isoforms/comparison_tss_5Pseq_fantom.csv",
            pdir=DIRECTORY)

rule consensus_sites:
    input:
        internal=expand("{{dir}}/mouse/5Pseq/tss/summary/{sample}_{replicate}.summary.counts.csv",
            replicate=REPLICATES,
            sample=SAMPLES),
        fantom=expand("{{dir}}/mouse/fantom/tss/summary/{sample}.summary.counts.csv",
            sample=FANTOM),
    output:
        collapsed="{dir}/combined/mouse_combined_consensus_sites.csv",
        summary="{dir}/combined/mouse_consensus_sites_summary.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_consensus.pl\
            --dir {wildcards.dir}/mouse/5Pseq/tss/summary,{wildcards.dir}/mouse/fantom/tss/summary \
            --collapsed {output.collapsed} \
            --summarised {output.summary} \
            --maxdist 12\
            --minbc 10\
            --suffix .summary.counts.csv \
        """

rule map_sample_sites:
    input:
        summary="{dir}/mouse/5Pseq/tss/summary/{sample}_{replicate}.summary.counts.csv",
        collapsed="{dir}/combined/mouse_combined_consensus_sites.csv",
    output:
        isoform="{dir}/mouse/5Pseq/tss/isoforms/{sample}_{replicate}.isoforms.csv",
        gene="{dir}/mouse/5Pseq/tss/genes/{sample}_{replicate}.genes.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_sample.pl\
            --insummary {input.summary} \
            --iniso {input.collapsed} \
            --outiso {output.isoform} \
            --outgene {output.gene} \
        """

rule map_fantom_sites:
    input:
        summary="{dir}/mouse/fantom/tss/summary/{fantom}.summary.counts.csv",
        collapsed="{dir}/combined/mouse_combined_consensus_sites.csv",
    output:
        isoform="{dir}/mouse/fantom/tss/isoforms/{fantom}.isoforms.csv",
        gene="{dir}/mouse/fantom/tss/genes/{fantom}.genes.csv",
    shell:
        """
        perl ~/analysis/tss/common-scripts/collapse_tss_sample.pl\
            --insummary {input.summary} \
            --iniso {input.collapsed} \
            --outiso {output.isoform} \
            --outgene {output.gene} \
        """

rule compare_sites:
    input:
        internal=expand("{{dir}}/mouse/5Pseq/tss/isoforms/{sample}_{replicate}.isoforms.csv",
            replicate=REPLICATES,
            sample=SAMPLES),
        fantom=expand("{{dir}}/mouse/fantom/tss/isoforms/{fantom}.isoforms.csv",
            fantom=FANTOM)
    output:
        comparison="{dir}/combined/tss/isoforms/comparison_tss_5Pseq_fantom.csv",
    shell:
        """
        Rscript ~/analysis/tss/mouse/aggregate_mESC_isoforms.r\
            --mousedir {wildcards.dir}/mouse/fantom/tss/isoforms \
            --fantomdir {wildcards.dir}/mouse/5Pseq/tss/isoforms \
            --odir {wildcards.dir}/combined/tss/isoforms \
            --suffix .isoforms.csv \
        """
