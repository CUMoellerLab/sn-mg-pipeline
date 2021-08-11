from os.path import splitext

rule index_contigs_bt2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.to_sample,
                                                binning_df)
    output:
        temp(multiext("output/binning/indexed/{to_sample}",
                      ".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2"))
    log:
        "output/logs/binning/bowtie2-build.{to_sample}.log"
    conda:
        "../env/bowtie2.yaml"
    params:
        bt2b_command = config['params']['bowtie2']['bt2b_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
        indexbase = "output/binning/indexed/{to_sample}"
    threads:
        config['threads']['bowtie2_build']
    shell:
        """
        {params.bt2b_command} --threads {threads} \
        {input.contigs} {params.indexbase} 2> {log} 1>&2
        """

rule map_reads_bt2:
    """
    Maps reads to contig files using bowtie2.
    """
    input:
        reads = lambda wildcards: expand("output/filtered/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.from_sample,
                                         read=['1','2']),
        db=rules.index_contigs_bt2.output
    output:
        aln=temp("output/binning/mapped_reads/{from_sample}.{to_sample}.bam")
    params:
        ref="output/binning/indexed/{to_sample}",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['map_reads']
    benchmark:
        "output/benchmarks/map_reads_bt2/{from_sample}.{to_sample}.benchmark.txt"
    log:
        "output/logs/map_reads_bt2/{from_sample}.{to_sample}.bowtie.log"
    shell:
        """
        # Map reads against reference genome
        {params.bt2_command} {params.extra} -p {threads} -x {params.ref} \
          -1 {input.reads[0]} -2 {input.reads[1]} \
          2> {log} | samtools view -bS - > {output.aln}

        """

rule sort_bam_bt2:
    """
    Sorts a bam file.
    """
    input:
        aln=rules.map_reads_bt2.output.aln
    output:
        bam="output/binning/bt2/mapped_reads/{from_sample}.{to_sample}.sorted.bam"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/sort_bam/bt2/{from_sample}.{to_sample}.sorted.txt"
    log:
        "output/logs/sort_bam/bt2/{from_sample}.{to_sample}.sort.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}

        """


rule index_contigs_minimap2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.to_sample,
                                                binning_df)
    output:
        index=temp("output/binning/indexed/{to_sample}.mmi")
    log:
        "output/logs/binning/minimap2.index.{to_sample}.log"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['minimap2_index']
    shell:
        """
        minimap2 -d {output.index} {input.contigs} -t {threads} 2> {log} 1>&2

        """

rule map_reads_minimap2:
    """
    Maps reads to contig files using minimap2.
    """
    input:
        reads = lambda wildcards: expand("output/filtered/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.from_sample,
                                         read=['1','2']),
        db=rules.index_contigs_minimap2.output.index
    output:
        aln=temp("output/binning/map_reads_minimap2/{from_sample}.{to_sample}.bam")
    params:
        k=config['params']['minimap2']['k']
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['minimap2_map_reads']
    benchmark:
        "output/benchmarks/map_reads_minimap2/{from_sample}.{to_sample}.benchmark.txt"
    log:
        "output/logs/map_reads_minimap2/{from_sample}.{to_sample}.bowtie.log"
    shell:
        """
        # Map reads against contigs
        minimap2 -a {input.db} {input.reads} -t {threads} -K {params.k} \
        2> {log} | samtools view -bS - > {output.aln}

        """

rule sort_bam_minimap2:
    """
    Sorts a bam file.
    """
    input:
        aln=rules.map_reads_minimap2.output.aln
    output:
        bam="output/binning/minimap2/mapped_reads/{from_sample}.{to_sample}.sorted.bam"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/sort_bam/minimap2/{from_sample}.{to_sample}.sorted.txt"
    log:
        "output/logs/sort_bam/minimap2/{from_sample}.{to_sample}.sort.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}

        """
