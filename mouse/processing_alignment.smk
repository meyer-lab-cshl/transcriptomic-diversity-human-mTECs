# how to run from ./:
# snakemake -s processing_alignment.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json \
# --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources} \
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['mESC1', 'mESC2', 'mESC3']
DIRECTORY='/sonas-hs/meyer/hpc/home/hmeyer/data/tss/mouse/5Pseq'

rule all:
    input:
        expand("{pdir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.combined.bedgraph",
            replicate=[1, 2],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/tss/combined/all_mESCs.positions.csv",
            pdir=DIRECTORY)

rule unique_reads:
    input:
        sam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.out.bam",
    output:
        unique="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.sam",
    shell:
        """
        samtools view -q 255 -h {input.sam} > {output.unique}
        """

rule strand_reads:
    input:
        sam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.sam",
    output:
        plus="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.plus.sam",
        minus="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.minus.sam",
    shell:
        """
        python processing_alignment/divide_by_strand.py \
            --samfile {input.sam} \
            --outprefix {wildcards.dir}/alignments/{wildcards.sample}_{wildcards.replicate}_Aligned.sortedByCoord.unique
        """

rule convert_bam:
    input:
        sam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.{strand}.sam",
    output:
        bam="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.{strand}.bam",
        sort="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.{strand}.sorted.bam",
        bai="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.{strand}.sorted.bam.bai",
    shell:
        """
        samtools view -h -b -S {input.sam} > {output.bam}
        samtools sort -m 3G {output.bam} {wildcards.dir}/alignments/{wildcards.sample}_{wildcards.replicate}_Aligned.sortedByCoord.unique.{wildcards.strand}.sorted
        samtools index {output.sort}
        """

rule make_bedgraph:
    input:
        sort="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.{strand}.sorted.bam",
        bai="{dir}/alignments/{sample}_{replicate}_Aligned.sortedByCoord.unique.{strand}.sorted.bam.bai",
    output:
        bedgraph="{dir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.{strand}.bedgraph",
    shell:
        """
        genomeCoverageBed -5 -dz -ibam {input.sort} > {output.bedgraph}
        """

rule merge_bedgraph:
    input:
        plus="{dir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.plus.bedgraph",
        minus="{dir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.minus.bedgraph",
    output:
        bedgraph="{dir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.combined.bedgraph",
        upload="{dir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.combined.ucsc.bedgraph",
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
        counts="{dir}/bedgraphs/{sample}_{replicate}_5prime_nz_readdepth.combined.bedgraph",
    output:
        bedgraph="{dir}/tss/bedgraphs/{sample}_{replicate}.bedgraph",
        counts="{dir}/tss/summary/{sample}_{replicate}.summary.counts.csv",
        positions="{dir}/tss/raw_positions/{sample}_{replicate}.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/process_counts.r \
            --type bedgraph \
            --species mouse \
            --ifile {input.counts} \
            --odir {wildcards.dir}/tss \
            --sample {wildcards.sample}_{wildcards.replicate}
        """

rule combine_counts:
    input:
        positions=expand("{{dir}}/tss/raw_positions/{sample}_{replicate}.positions.csv",
            sample=SAMPLE,
            replicate=[1,2]),
    output:
        positions="{dir}/tss/combined/all_mESCs.positions.csv",
    shell:
        """
        Rscript ~/analysis/tss/common-scripts/combine_counts.r \
            --indir {wildcards.dir}/tss/raw_positions \
            --ofile {output.positions} \
            --suffix .positions.csv \
            --verbose
        """
