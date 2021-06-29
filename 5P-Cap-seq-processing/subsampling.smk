# how to run from ./:
# snakemake -s subsampling.smk --profile uge
#
# --jobs 500 --latency-wait 90 --cluster-config config/cluster.json --cluster 'qsub {cluster.nodes} -N {cluster.name}  -l {cluster.resources}
# -o {cluster.output} -e {cluster.error} -cwd' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_subsample.yml"

checkpoint generate_subsampling_matrix:
    input:
        expand("{dir}/subsampling/overview.dedup.unique.reads.txt",
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

def subsample_reads(wildcards):
    filename = checkpoints.generate_subsampling_matrix.get().output[0]
    directory = config['directory']
    files = []
    with open(filename) as f:
        lines = f.readlines()
        for l in lines:
            columns = l.strip().split("\t")
            sample = columns[0]
            reads = columns[1]
            print(f"{directory} {reads} {sample}")
            files.append(f"{directory}/subsampling/reads{reads}/Paraclu_BED/{sample}_TPM_paraclu_simplified_bed.txt")
    return files

rule all:
    input:
        expand("{pdir}/deduplicated/overview.dedup.unique.reads.txt",
            pdir=config['directory']),
        #expand("{pdir}/subsampling/reads{reads}/Paraclu_BED/{sample}_TPM_paraclu_simplified_bed.txt",
        #expand("{pdir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam",
        #    pdir=config['directory'],
        #    reads=[100000, 500000, 2000000],
        #    sample=config['samples']),
        subsample_reads,
        #expand("{pdir}/subsampling/reads{reads}/Merged_CTSS_BED/All_thymus_merged.txt",
        #    pdir=config['directory'],
        #    reads=[100000, 500000, 2000000])

rule count_reads:
    input:
        sort="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam",
        index="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam.bai",
    output:
        "{dir}/subsampling/{sample}_Aligned.sortedByCoord.dedup.unique.reads.txt"
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
        expand("{{dir}}/subsampling/{sample}_Aligned.sortedByCoord.dedup.unique.reads.txt",
            sample=config['samples'])
    output:
        "{dir}/subsampling/overview.dedup.unique.reads.txt"
    conda:
        "envs/subsampling.yaml"
    shell:
        """
        cat {input} >> {output}
        """

rule subsample:
    input:
        counts="{dir}/subsampling/subsampling_matrix.txt",
        #counts="{dir}/subsampling/overview.dedup.unique.reads.txt",
        #bam="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam",
        #bai="{dir}/deduplicated/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam.bai",
    output:
        "{dir}/subsampling/reads{reads}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam",
        #directory("{dir}/subsampling/reads{reads}")
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
            -r {wildcards.reads} \
            -t {threads} \
            -o {output}
        """

rule ctss_normalize:
    input:
        expand("{{dir}}/subsampling/reads{{reads}}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam",
            sample=config['samples'])
        #glob_wildcards("{{dir}}/subsampling/reads{{reads}}/{sample}_Aligned.sortedByCoord.dedup.unique.fwd.bam")
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
    params:
        samples = config['samples']
    script:
        "subsampling/ctss_normalize.R"

rule ctss_format:
    input:
        tpm="{dir}/subsampling/reads{reads}/ctss/TPM.csv",
        tagcounts="{dir}/subsampling/reads{reads}/ctss/tagCount.csv",
    output:
        expand("{{dir}}/subsampling/reads{{reads}}/for_Paraclu/{sample}_TPM_for_Paraclu.txt",
            sample=config['samples']),
        expand("{{dir}}/subsampling/reads{{reads}}/CAGEr_out/DFs/{sample}_Counts_and_TPM.txt",
            sample=config['samples'])
    conda:
        "envs/ctss.yaml"
    log:
        "{dir}/subsampling/reads{reads}/log/ctss_format.out",
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



