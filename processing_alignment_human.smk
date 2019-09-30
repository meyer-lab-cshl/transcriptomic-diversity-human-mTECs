# how to run from ./:
# snakemake -s processing_alignment_human.smk --jobs 5000 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_processing_alignment.yml"

SAMPLE=['pt87-hi', 'pt87-lo', 'pt212-hi', 'pt212-lo', 'pt214-hi', 'pt214-lo',
        'pt221-hi', 'pt221-lo', 'pt226-hi', 'pt226-lo']
DIRECTORY='/mnt/grid/meyer/hpc/home/data/common/tss'

rule all:
    input:
        expand("{pdir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.{read}.sam",
            read=['fwd', 'rev'],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/3_tss_data/raw_positions/single_samples/{sample}_fwd.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/3_tss_data/raw_positions/all_mTECs_fwd.positions.csv",
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.{strand}.bam",
            strand=['plus', 'minus'],
            pdir=DIRECTORY,
            sample=SAMPLE),
        expand("{pdir}/2_alignments/{sample}_5prime_nz_readdepth.combined.bedgraph",
            pdir=DIRECTORY,
            sample=SAMPLE),

rule unique_reads:
    input:
        sam="{dir}/2_alignments/{sample}_Aligned.out.sam",
    output:
        unique="{dir}/2_alignments/{sample}_Aligned.out.unique.sam",
    shell:
        """
        samtools view -q 255 -S -h {input.sam} > {output.unique}
        """

rule add_feature:
    input:
        unique="{dir}/2_alignments/{sample}_Aligned.out.unique.sam",
        info="{dir}/1_raw_data/5Pseq_processed/molbc/{sample}.molbc.tsv",
    output:
        mbc="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.sam",
    shell:
        """
        python processing_alignment/add_feature.py \
            --sam {input.unique} \
            --info {input.info} \
            --infoname molecular_bc \
            --flagname M \
            --out {output.mbc}
        """

rule split_reads:
    input:
        mbc="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.sam",
    output:
        fwd="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.fwd.sam",
        rev="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.rev.sam",
    shell:
        """
        samtools view -h -f 0x40 -S {input.mbc} > {output.fwd}
        samtools view -h -f 0x80 -S {input.mbc} > {output.rev}
        """

rule collapse_reads:
    input:
        sam="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.{read}.sam",
    output:
        collapsed="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.{read}.sam",
        mc_count="{dir}/2_alignments/collapsing_{read}/{sample}_collapsed_countsPerMol.tsv",
        pos_count="{dir}/2_alignments/collapsing_{read}/{sample}_collapsed_countsPerPos.tsv",
        info="{dir}/2_alignments/collapsing_{read}/{sample}_collapsed.info",
    shell:
        """
        python processing_alignment/collapse_reads.py \
            --samfile {input.sam} \
            --out-collapsed {output.collapsed} \
            --out-molecule-count {output.mc_count} \
            --out-position-count {output.pos_count} \
            --out-info {output.info}
        """

rule process_counts:
    input:
        counts="{dir}/2_alignments/collapsing_fwd/{sample}_collapsed_countsPerMol.tsv",
    output:
        bedgraph="{dir}/3_tss_data/bedgraphs/single_samples/{sample}_fwd.bedgraph",
        counts="{dir}/3_tss_data/summary/single_samples/{sample}_fwd.summary.counts.csv",
        positions="{dir}/3_tss_data/raw_positions/single_samples/{sample}_fwd.positions.csv",
    shell:
        """
        Rscript processing_alignment/process_counts.r \
            --ifile {input.counts} \
            --odir {wildcards.dir}/3_tss_data \
            --sample {wildcards.sample}_fwd
        """

rule combine_counts:
    input:
        positions=expand("{{dir}}/3_tss_data/raw_positions/single_samples/{sample}_fwd.positions.csv",
            sample=SAMPLE),
    output:
        positions="{dir}/3_tss_data/raw_positions/all_mTECs_fwd.positions.csv",
    shell:
        """
        Rscript processing_alignment/combine_counts.r \
            --directory {wildcards.dir}/3_tss_data/raw_positions \
            --suffix _fwd.positions.csv \
            --ofile {output.positions} \
            --verbose
        """

rule strand_reads:
    input:
        sam="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.sam",
    output:
        plus="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.plus.sam",
        minus="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.minus.sam",
    shell:
        """
        python processing_alignment/divide_by_strand.py \
            --samfile {input.sam} \
            --outprefix {wildcards.dir}/2_alignments/{wildcards.sample}_Aligned.out.unique.mbc.collapsed.fwd
        """

rule convert_bam:
    input:
        sam="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.{strand}.sam",
    output:
        bam="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.{strand}.bam",
        sort="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.{strand}.sorted.bam",
    shell:
        """
        samtools view -h -b -S {input.sam} > {output.bam}
        samtools sort -f {output.bam} {output.sort}
        samtools index {output.sort}
        """

rule make_bedgraph:
    input:
        sort="{dir}/2_alignments/{sample}_Aligned.out.unique.mbc.collapsed.fwd.{strand}.sorted.bam",
    output:
        bedgraph="{dir}/2_alignments/{sample}_5prime_nz_readdepth.{strand}.bedgraph",
    shell:
        """
        genomeCoverageBed -5 -dz -ibam {input.sort} > {output.bedgraph}
        """

rule merge_bedgraph:
    input:
        plus="{dir}/2_alignments/{sample}_5prime_nz_readdepth.plus.bedgraph",
        minus="{dir}/2_alignments/{sample}_5prime_nz_readdepth.minus.bedgraph",
    output:
        bedgraph="{dir}/2_alignments/{sample}_5prime_nz_readdepth.combined.bedgraph",
        upload="{dir}/2_alignments/{sample}_5prime_nz_readdepth.combined.ucsc.bedgraph",
    shell:
        """
        python processing_alignment/merge_bedgraph.py \
            --plus {input.plus} \
            --minus {input.minus} \
            --outfile {output.bedgraph}
        echo 'track type=bedGraph name="" description="BedGraph format" visibility=full color=200,100,0' > {output.upload}
        cat {output.bedgraph} >> {output.upload}
        """

