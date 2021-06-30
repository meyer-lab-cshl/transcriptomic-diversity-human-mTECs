# how to run from ./:
# snakemake -s subsampling.smk --profile uge --use-conda --until all

configfile: "config/config_subsample.yml"

checkpoint generate_subsampling_matrix:
    input:
        expand("{dir}/subsampling/overview.reads.txt",
            dir=config['directory'])
    output:
        expand("{dir}/subsampling/subsampling_matrix.txt",
            dir=config['directory'])
    log:
        expand("{dir}/subsampling/log/subsampling_matrix.out",
            dir=config['directory'])
    conda:
        "envs/subsampling.yaml"
    script:
        "subsampling/subsampling_matrix.R"


def subsample_reads(wilcards):
    filename = checkpoints.generate_subsampling_matrix.get().output[0]
    directory = config['directory']
    files = []
    with open(filename) as f:
        lines = f.readlines()
        for l in lines[1:]:
            columns = l.strip().split("\t")
            sample = columns[0]
            reads = columns[1]
            files.append(f"{directory}/subsampling/reads{reads}/Paraclu_BED/{sample}_TPM_paraclu_simplified_bed.txt")
    return files

def ctss_files(wildcards):
    filename = checkpoints.generate_subsampling_matrix.get().output[0]
    directory = config['directory']
    files = []
    with open(filename) as f:
        lines = f.readlines()
        for l in lines[1:]:
            columns = l.strip().split("\t")
            sample = columns[0]
            reads = columns[1]
            if reads == wildcards.reads:
                files.append(f"{directory}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam")
    return files

rule all:
    input:
        expand("{pdir}/subsampling/overview.reads.txt",
            pdir=config['directory']),
        expand("{dir}/subsampling/subsampling_matrix.txt",
            dir=config['directory']),
        subsample_reads,
        expand("{pdir}/subsampling/subsampling_cluster_summary.txt",
            pdir=config['directory']),


rule count_reads:
    input:
        sort="{dir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        index="{dir}/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "{dir}/subsampling/{sample}_Aligned.sortedByCoord.reads.txt"
    conda:
        "envs/subsampling.yaml"
    shell:
        """
        samtools idxstats {input.sort} | \
            cut -f3 | \
            awk -v s={wildcards.sample} 'BEGIN {{total=0}} {{total += $1}} END {{print s,"\t",total}}' > {output}
        """


rule overview_count_reads:
    input:
        expand("{{dir}}/subsampling/{sample}_Aligned.sortedByCoord.reads.txt",
            sample=config['samples'])
    output:
        "{dir}/subsampling/overview.reads.txt"
    conda:
        "envs/subsampling.yaml"
    shell:
        """
        cat {input} >> {output}
        """

rule subsample:
    input:
        counts="{dir}/subsampling/subsampling_matrix.txt",
        bam="{dir}/alignments/{sample}_Aligned.sortedByCoord.out.bam",
        bai="{dir}/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.bam",
    wildcard_constraints:
        reads="\d+"
    threads: 2
    resources:
        mem_mb = 1000
    conda:
        "envs/subsampling.yaml"
    shell:
        """
        subsampling/subsample.sh \
            -c {input.counts} \
            -s {wildcards.sample} \
            -r {wildcards.reads} \
            -i {input.bam} \
            -t {threads} \
            -o {output}
        """

rule dedup:
    input:
        bam="{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.bam",
        index="{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.bam.bai",
    output:
        bam="{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.bam",
    resources:
        mem_mb = 12000
    threads: 4
    conda:
        "envs/subsampling.yaml"
    shell:
        """
        umi_tools dedup \
            -I {input.bam} \
            -L {wildcards.dir}/deduplicated/extract.log \
            -S {output.bam} \
            --paired \
            --output-stats={wildcards.dir}/subsampling/reads{wildcards.reads}/{wildcards.sample}_Aligned.sortedByCoord.stats_
        """

rule process_dedup:
    input:
        bam="{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.bam",
    output:
        fwd=temp("{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.sam"),
        unique=temp("{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.bam"),
        sort="{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam",
        index="{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam.bai",
    resources:
        mem_mb = 12000
    threads: 4
    conda:
        "envs/subsampling.yaml"
    shell:
        """
        samtools view -b -q 255 {input.bam} > {output.unique}
        samtools view -h -f 0x40 {output.unique} > {output.fwd}
        samtools sort {output.fwd} -o {output.sort} -@ {threads} -m 20G
        samtools index {output.sort} -@ {threads}
        """


rule ctss_normalize:
    input:
        ctss_files
    output:
        powerlaw="{dir}/subsampling/reads{reads}/ctss/PowerLaw.pdf",
        powerlawnorm="{dir}/subsampling/reads{reads}/ctss/PowerLaw_Normalized.pdf",
        tpm="{dir}/subsampling/reads{reads}/ctss/TPM.csv",
        tagcounts="{dir}/subsampling/reads{reads}/ctss/tagCount.csv",
    wildcard_constraints:
        reads="\d+"
    conda:
        "envs/ctss.yaml"
    log:
        "{dir}/subsampling/reads{reads}/log/ctss_normalize.out",
    script:
        "subsampling/ctss_normalize.R"


rule ctss_format:
    input:
        tpm="{dir}/subsampling/reads{reads}/ctss/TPM.csv",
        tagcounts="{dir}/subsampling/reads{reads}/ctss/tagCount.csv",
    output:
        "{dir}/subsampling/reads{reads}/for_Paraclu/{sample}_TPM_for_Paraclu.txt",
        "{dir}/subsampling/reads{reads}/CAGEr_out/DFs/{sample}_Counts_and_TPM.txt",
    conda:
        "envs/ctss.yaml"
    #log:
    #    "{dir}/subsampling/reads{reads}/log/ctss_format.out",
    script:
        "subsampling/ctss_format.R"


rule call_clusters:
    input:
        "{dir}/subsampling/reads{reads}/for_Paraclu/{sample}_TPM_for_Paraclu.txt"
    output:
        "{dir}/subsampling/reads{reads}/Paraclu_Cluster/{sample}_TPM_paraclu.txt"
    params:
        min_TPM=2,
    threads: 2
    resources:
        mem_mb = 1000
    conda:
        "envs/clusters.yaml"
    shell:
        """
        #Paraclu - tag clusters
        paraclu {params.min_TPM} {input} > {output}
        """

rule trim_clusters:
    input:
        "{dir}/subsampling/reads{reads}/Paraclu_Cluster/{sample}_TPM_paraclu.txt"
    output:
        "{dir}/subsampling/reads{reads}/Paraclu_Trim/{sample}_TPM_paraclu_simplified.txt"
    params:
        max_length=20
    threads: 2
    resources:
        mem_mb = 1000
    conda:
        "envs/clusters.yaml"
    shell:
        """
        #Filter Paraclu - apply max length, min expression thresholds
        paraclu-cut -l {params.max_length} {input} > {output}
        """

rule process_clusters:
    input:
        "{dir}/subsampling/reads{reads}/Paraclu_Trim/{sample}_TPM_paraclu_simplified.txt"
    output:
        "{dir}/subsampling/reads{reads}/Paraclu_BED/{sample}_TPM_paraclu_simplified_bed.txt"
    resources:
        mem_mb = 1000
    shell:
        """
        #Process paraclu output, export as BED format
        awk '{{print $1 "\\t" $3 "\\t" $4 "\\t" $2}}' {input} > {output}_tmp1
        sed -i "s/$/\\t{wildcards.sample}\\t{wildcards.sample}_/" {output}_tmp1
        awk '{{ print $0, NR }}' {output}_tmp1 > {output}_tmp2
        awk '{{new_var=$6$7; print $0, new_var}}' {output}_tmp2  > {output}_tmp3
        awk '{{print $1 "\\t" $2 "\\t" $3 "\\t" $5 "\\t" $8 "\\t" $4}}' {output}_tmp3  \
            > {output}
        rm {output}_tmp*
        """

rule consensus_clusters:
    input:
        expand("{{dir}}/subsampling/reads{{reads}}/Paraclu_BED/{sample}_TPM_paraclu_simplified_bed.txt",
            sample=config['samples'])
    output:
        "{dir}/subsampling/reads{reads}/Merged_CTSS_BED/All_thymus_merged.txt"
    conda:
        "envs/clusters.yaml"
    shell:
        """
        bedops --everything {input} > {output}_tmp_union_bed
        sort-bed {output}_tmp_union_bed > {output}_tmp_sorted_union_bed
        bedtools merge -i {output}_tmp_sorted_union_bed -s -d 20 -c 4,5,6 \
                  -o distinct,distinct,distinct > {output}
        rm {output}_tmp*
        """

rule cluster_summary:
    input:
        subsample_reads
    output:
        "{dir}/subsampling/subsampling_cluster_summary.txt"
    shell:
        """
        echo "sample\treads\tclusters" > {output}
        for file in {input}; do
            filename=${{file##*/}}
            sample=${{filename%%_*}}
            tmp=${{file##*/reads}}
            reads=${{tmp%%/*}}
            clusters=$(wc -l $file | cut -d " "  -f 1)
            echo "$sample\t$reads\t$clusters" >> {output}
        done
        """

