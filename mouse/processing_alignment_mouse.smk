# how to run from ./:
# snakemake -s processing_alignment.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['mESC1', 'mESC2', 'mESC3']
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss'

rule all:
    input:
        expand("{pdir}/mouse/2_alignments/{sample}_{replicate}_Aligned.out.unique.sam",
            replicate=[1, 2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/mouse/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.bam",
            replicate=[1, 2],
            strand=['plus', 'minus'],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/mouse/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.{strand}.bedgraph",
            strand=['plus', 'minus'],
            replicate=[1, 2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/mouse/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.combined.bedgraph",
            replicate=[1, 2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/mouse/3_tss_data/raw_positions/{sample}_{replicate}.positions.csv",
            replicate=[1, 2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/combined/3_tss_data/raw_positions/all_mESCs.positions.csv",
            pdir=DIRECTORY)

rule unique_reads:
    input:
        sam="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.sam",
    output:
        unique="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.sam",
    shell:
        """
        samtools view -q 255 -S -h {input.sam} > {output.unique}
        """

rule strand_reads:
    input:
        sam="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.sam",
    output:
        plus="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.plus.sam",
        minus="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.minus.sam",
    shell:
        """
        python processing_alignment/divide_by_strand.py \
            --samfile {input.sam} \
            --outprefix {wildcards.dir}/2_alignments/{wildcards.sample}_{wildcards.replicate}_Aligned.out.unique
        """

rule convert_bam:
    input:
        sam="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.sam",
    output:
        bam="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.bam",
        sort="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.sorted.bam",
        bai="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.sorted.bam.bai",
    shell:
        """
        samtools view -h -b -S {input.sam} > {output.bam}
        samtools sort -m 3G {output.bam} {wildcards.dir}/2_alignments/{wildcards.sample}_{wildcards.replicate}_Aligned.out.unique.{wildcards.strand}.sorted
        samtools index {output.sort}
        """

rule make_bedgraph:
    input:
        sort="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.sorted.bam",
        bai="{dir}/2_alignments/{sample}_{replicate}_Aligned.out.unique.{strand}.sorted.bam.bai",
    output:
        bedgraph="{dir}/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.{strand}.bedgraph",
    shell:
        """
        genomeCoverageBed -5 -dz -ibam {input.sort} > {output.bedgraph}
        """

rule merge_bedgraph:
    input:
        plus="{dir}/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.plus.bedgraph",
        minus="{dir}/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.minus.bedgraph",
    output:
        bedgraph="{dir}/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.combined.bedgraph",
        upload="{dir}/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.combined.ucsc.bedgraph",
    shell:
        """
        python processing_alignment/merge_bedgraph.py \
            --plus {input.plus} \
            --minus {input.minus} \
            --outfile {output.bedgraph}
        echo 'track type=bedGraph name="" description="BedGraph format" visibility=full color=200,100,0' > {output.upload}
        cat {output.bedgraph} >> {output.upload}
        """

rule process_counts:
    input:
        counts="{dir}/2_alignments/{sample}_{replicate}_5prime_nz_readdepth.combined.bedgraph",
    output:
        bedgraph="{dir}/3_tss_data/bedgraphs/{sample}_{replicate}.bedgraph",
        counts="{dir}/3_tss_data/summary/{sample}_{replicate}.summary.counts.csv",
        positions="{dir}/3_tss_data/raw_positions/{sample}_{replicate}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --type bedgraph \
            --species mouse \
            --ifile {input.counts} \
            --odir {wildcards.dir}/3_tss_data \
            --sample {wildcards.sample}_{wildcards.replicate}
        """

rule combine_counts:
    input:
        positions=expand("{{dir}}/mouse/3_tss_data/raw_positions/{sample}_{replicate}.positions.csv",
            sample=SAMPLE,
            replicate=[1,2]),
    output:
        positions="{dir}/combined/3_tss_data/raw_positions/all_mESCs.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/mouse/3_tss_data/raw_positions \
            --odir {wildcards.dir}/combined/3_tss_data/raw_positions \
            --suffix .positions.csv \
            --verbose
        """
