# how to run from ./:
# snakemake -s processing_alignment.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
        'pt221-hi', 'pt221-lo', 'pt226-hi', 'pt226-lo']
DIRECTORY='/mnt/grid/meyer/hpc/home/data/hmeyer/data/tss/human/5Pseq'

rule all:
    input:
        expand("{pdir}/collapsed/{sample}_collapsed_{reads}.info",
            reads=['fwd', 'rev'],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.{reads}.sam",
            reads=['fwd', 'rev'],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/bedgraphs/{sample}_5prime_nz_readdepth.combined.bedgraph",
            pdir=DIRECTORY,
            sample=SAMPLE)

rule unique_reads:
    input:
        bam="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.bam",
    output:
        unique="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.bam",
    shell:
        """
        samtools view -b -q 255 {input.bam} > {output.unique}
        """

rule split_reads:
    input:
        unique="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.bam",
    output:
        fwd="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.sam",
        rev="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.rev.sam",
    shell:
        """
        samtools view -h -f 0x40 {input.unique} > {output.fwd}
        samtools view -h -f 0x80 {input.unique} > {output.rev}
        """

rule collapse_reads:
    input:
        sam="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.{read}.sam",
    output:
        collapsed="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.{read}.sam",
        pos_count="{dir}/collapsed/{sample}_collapsed_{read}_countsPerPos.tsv",
        info="{dir}/collapsed/{sample}_collapsed_{read}.info",
    shell:
        """
        python processing_alignment/collapse_reads.py \
            --samfile {input.sam} \
            --out-collapsed {output.collapsed} \
            --out-position-count {output.pos_count} \
            --out-info {output.info}
        """

rule strand_reads:
    input:
        sam="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.sam",
    output:
        plus="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.plus.sam",
        minus="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.minus.sam",
    shell:
        """
        python processing_alignment/divide_by_strand.py \
            --samfile {input.sam} \
            --outprefix {wildcards.dir}/collapsed/{wildcards.sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd
        """

rule convert_bam:
    input:
        sam="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.{strand}.sam",
    output:
        bam="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.{strand}.bam",
        sort="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.{strand}.sorted.bam",
    shell:
        """
        samtools view -h -b -S {input.sam} > {output.bam}
        samtools sort -f {output.bam} {output.sort}
        samtools index {output.sort}
        """

rule make_bedgraph:
    input:
        sort="{dir}/collapsed/{sample}_Aligned.sortedByCoord.dedup.unique.collapsed.fwd.{strand}.sorted.bam",
    output:
        bedgraph="{dir}/bedgraphs/{sample}_5prime_nz_readdepth.{strand}.bedgraph",
    shell:
        """
        genomeCoverageBed -5 -dz -ibam {input.sort} > {output.bedgraph}
        """

rule merge_bedgraph:
    input:
        plus="{dir}/bedgraphs/{sample}_5prime_nz_readdepth.plus.bedgraph",
        minus="{dir}/bedgraphs/{sample}_5prime_nz_readdepth.minus.bedgraph",
    output:
        bedgraph="{dir}/bedgraphs/{sample}_5prime_nz_readdepth.combined.bedgraph",
        upload="{dir}/bedgraphs/{sample}_5prime_nz_readdepth.combined.ucsc.bedgraph",
    shell:
        """
        python processing_alignment/merge_bedgraph.py \
            --spikeins M24537,X04603,X17013 \
            --plus {input.plus} \
            --minus {input.minus} \
            --outfile {output.bedgraph}
        echo 'track type=bedGraph name="" description="BedGraph format" visibility=full color=200,100,0' > {output.upload}
        cat {output.bedgraph} >> {output.upload}
        """

