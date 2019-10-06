# how to run from ./:
# snakemake -s processing_fantom.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

#configfile: "config/config_processing_alignment.yml"

DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss'
SAMPLE_MOUSE=['ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1.mm9.ctss',
'ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep2.CNhs14109.14362-155I6.mm9.ctss']
SAMPLE_HUMAN, = glob_wildcards("/sonas-hs/meyer/hpc/home/hmeyer/data/tss/fantom/human/{sample}.bed.gz")

rule all:
    input:
        expand("{pdir}/fantom/mouse/3_tss_data/raw_positions/{sample}.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE_MOUSE),
        expand("{pdir}/combined/3_tss_data/raw_positions/all_mESC_46C.mm9.ctss.positions.csv",
            pdir=DIRECTORY),
        expand("{pdir}/fantom/human/3_tss_data/raw_positions/{sample}.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE_HUMAN),
        expand("{pdir}/combined/3_tss_data/raw_positions/all_tissues.hg19.ctss.positions.csv",
            pdir=DIRECTORY),

rule process_counts_mouse:
    input:
        counts="{dir}/fantom/mouse/mouse.timecourse.hCAGE/{sample}.bed.gz",
    output:
        bedgraph="{dir}/fantom/mouse/3_tss_data/bedgraphs/{sample}.bedgraph",
        counts="{dir}/fantom/mouse/3_tss_data/summary/{sample}.summary.counts.csv",
        positions="{dir}/fantom/mouse/3_tss_data/raw_positions/{sample}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --species mouse \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/fantom/mouse/3_tss_data \
            --sample {wildcards.sample}
        """

rule combine_counts_mouse:
    input:
        positions=expand("{{dir}}/fantom/mouse/3_tss_data/raw_positions/{sample}.positions.csv",
            sample=SAMPLE_MOUSE,
            replicate=[1,2]),
    output:
        positions="{dir}/combined/3_tss_data/raw_positions/all_mESC_46C.mm9.ctss.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/fantom/mouse/3_tss_data/raw_positions \
            --ofile {output.positions} \
            --suffix .mm9.ctss.positions.csv \
            --verbose
        """

rule process_counts_human:
    input:
        counts="{dir}/fantom/human/{sample}.bed.gz",
    output:
        bedgraph="{dir}/fantom/human/3_tss_data/bedgraphs/{sample}.bedgraph",
        counts="{dir}/fantom/human/3_tss_data/summary/{sample}.summary.counts.csv",
        positions="{dir}/fantom/human/3_tss_data/raw_positions/{sample}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --species human \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/fantom/human/3_tss_data \
            --sample {wildcards.sample}
        """

rule combine_counts_human:
    input:
        positions=expand("{{dir}}/fantom/human/3_tss_data/raw_positions/{sample}.positions.csv",
            sample=SAMPLE_HUMAN)
    output:
        positions="{dir}/combined/3_tss_data/raw_positions/all_tissues.hg19.ctss.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/fantom/human/3_tss_data/raw_positions \
            --odir {wildcards.dir}/combined/3_tss_data/raw_positions \
            --suffix .hg19.ctss.positions.csv \
            --fantom \
            --verbose
        """
