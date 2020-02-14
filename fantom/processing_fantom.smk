# how to run from ./:
# snakemake -s processing_fantom.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

#configfile: "config/config_processing_alignment.yml"

DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss'
SAMPLE_MOUSE=['ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1.mm9.ctss',
'ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep2.CNhs14109.14362-155I6.mm9.ctss']
SAMPLE_HUMAN, = glob_wildcards("/sonas-hs/meyer/hpc/home/hmeyer/data/tss/human/fantom/bed/{sample}.bed.gz")

rule all:
    input:
        expand("{pdir}/mouse/fantom/tss/raw_positions/{sample}.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE_MOUSE),
        expand("{pdir}/mouse/fantom/tss/combined/all_mESC_46C.mm9.ctss.positions.csv",
            pdir=DIRECTORY),
        expand("{pdir}/human/fantom/tss/raw_positions/{sample}.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE_HUMAN),
        expand("{pdir}/human/fantom/tss/combined/all_tissues.hg19.ctss.positions.csv",
            pdir=DIRECTORY),

rule process_counts_mouse:
    input:
        counts="{dir}/mouse/fantom/bed/mouse.timecourse.hCAGE/{sample}.bed.gz",
    output:
        bedgraph="{dir}/mouse/fantom/tss/bedgraphs/{sample}.bedgraph",
        counts="{dir}/mouse/fantom/tss/summary/{sample}.summary.counts.csv",
        positions="{dir}/mouse/fantom/tss/raw_positions/{sample}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --species mouse \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/mouse/fantom/tss \
            --sample {wildcards.sample}
        """

rule combine_counts_mouse:
    input:
        positions=expand("{{dir}}/mouse/fantom/tss/raw_positions/{sample}.positions.csv",
            sample=SAMPLE_MOUSE,
            replicate=[1,2]),
    output:
        positions="{dir}/mouse/fantom/tss/combined/all_mESC_46C.mm9.ctss.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/mouse/fantom/tss/raw_positions \
            --ofile {output.positions} \
            --suffix .mm9.ctss.positions.csv \
            --verbose
        """

rule process_counts_human:
    input:
        counts="{dir}/human/fantom/bed/{sample}.bed.gz",
    output:
        bedgraph="{dir}/human/fantom/tss/bedgraphs/{sample}.bedgraph",
        counts="{dir}/human/fantom/tss/summary/{sample}.summary.counts.csv",
        positions="{dir}/human/fantom/tss/raw_positions/{sample}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --species human \
            --type bed \
            --ifile {input.counts} \
            --odir {wildcards.dir}/human/fantom/tss \
            --sample {wildcards.sample}
        """

rule combine_counts_human:
    input:
        positions=expand("{{dir}}/human/fantom/tss/raw_positions/{sample}.positions.csv",
            sample=SAMPLE_HUMAN)
    output:
        positions="{dir}/human/fantom/tss/combined/all_tissues.hg19.ctss.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/human/fantom/tss/raw_positions \
            --odir {wildcards.dir}/human/fantom/tss/raw_positions \
            --suffix .hg19.ctss.positions.csv \
            --fantom \
            --verbose
        """
