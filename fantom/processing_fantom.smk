# how to run from ./:
# snakemake -p -s processing_fantom.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json
# --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

#configfile: "config/config_processing_alignment.yml"

DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss'
ANNOTATION='/sonas-hs/meyer/hpc/home/hmeyer/data/common/public/annotations/genome/liftover'
SAMPLE_MOUSE=['ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep1.CNhs14104.14357-155I1',
'ES-46C_embryonic_stem_cells_neuronal_differentiation_day00_biol_rep2.CNhs14109.14362-155I6']
SAMPLE_HUMAN, = glob_wildcards("/sonas-hs/meyer/hpc/home/hmeyer/data/tss/human/fantom/bed/GRCh37/{sample}.hg19.ctss.bed.gz")

rule all:
    input:
        expand("{pdir}/mouse/fantom/tss/combined/all_mESC_46C.positions.csv",
            pdir=DIRECTORY),
        #expand("{pdir}/human/fantom/tss/combined/all_tissues.positions.csv",
        #    pdir=DIRECTORY),

rule liftover_counts_mouse:
    input:
        positions="{dir}/mouse/fantom/bed/GRCm37/{sample}.mm9.ctss.bed.gz",
        chain=expand("{annodir}/mm9ToMm10.over.chain",
            annodir=ANNOTATION),
    output:
        liftover="{dir}/mouse/fantom/bed/GRCm38/{sample}.bed",
        unmapped="{dir}/mouse/fantom/bed/GRCm38/{sample}.unmapped.bed",
    shell:
        """
        liftOver {input.positions} {input.chain} \
            {output.liftover}.tmp {output.unmapped}.tmp
        awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2".."$3","$6,$5,$6}}' {output.liftover}.tmp \
            > {output.liftover}
        awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2".."$3","$6,$5,$6}}' {output.unmapped}.tmp \
            > {output.unmapped}
        rm {output.liftover}.tmp {output.unmapped}.tmp
        """

rule process_counts_mouse:
    input:
        counts="{dir}/mouse/fantom/bed/GRCm38/{sample}.bed",
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
        positions="{dir}/mouse/fantom/tss/combined/all_mESC_46C.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/mouse/fantom/tss/raw_positions \
            --ofile {output.positions} \
            --suffix .positions.csv \
            --verbose
        """


rule liftover_counts_human:
    input:
        positions="{dir}/human/fantom/bed/GRCh37/{sample}.hg19.ctss.bed.gz",
        chain=expand("{annodir}/hg19ToHg38.over.chain",
            annodir=ANNOTATION),
    output:
        liftover="{dir}/human/fantom/bed/GRCh38/{sample}.bed",
        unmapped="{dir}/human/fantom/bed/GRCh38/{sample}.unmapped.bed",
    shell:
        """
        liftOver {input.positions} {input.chain} \
            {output.liftover}.tmp {output.unmapped}.tmp
        awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2".."$3","$6,$5,$6}}' {output.liftover}.tmp \
            > {output.liftover}
        awk -v OFS="\\t" '{{print $1,$2,$3,$1":"$2".."$3","$6,$5,$6}}' {output.unmapped}.tmp \
            > {output.unmapped}
        rm {output.liftover}.tmp {output.unmapped}.tmp
        """

rule process_counts_human:
    input:
        counts="{dir}/human/fantom/bed/GRCh38/{sample}.bed",
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
        positions="{dir}/human/fantom/tss/combined/all_tissues.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/human/fantom/tss/raw_positions \
            --odir {wildcards.dir}/human/fantom/tss/combined \
            --suffix .positions.csv \
            --fantom \
            --verbose
        """


