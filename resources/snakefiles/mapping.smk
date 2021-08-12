from os.path import splitext

rule index_contigs_bt2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.to_sample,
                                                binning_df)
    output:
        temp(multiext("output/binning/bowtie2/indexed_contigs/{to_sample}",
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
        indexbase = "output/binning/bowtie2/indexed_contigs/{to_sample}"
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
        aln=temp("output/binning/bowtie2/mapped_reads/{from_sample}.MappedTo.{to_sample}.bam")
    params:
        ref="output/binning/bowtie2/indexed_contigs/{to_sample}",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['map_reads']
    benchmark:
        "output/benchmarks/map_reads_bt2/{from_sample}.MappedTo.{to_sample}.benchmark.txt"
    log:
        "output/logs/map_reads_bt2/{from_sample}.MappedTo.{to_sample}.bowtie.log"
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
        bam="output/binning/bowtie2/mapped_reads/{from_sample}.MappedTo.{to_sample}.sorted.bam"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/sort_bam/bowtie2/{from_sample}.MappedTo.{to_sample}.sorted.txt"
    log:
        "output/logs/sort_bam/bowtie2/{from_sample}.MappedTo.{to_sample}.sort.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}

        """


rule index_contigs_minimap2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.to_sample,
                                                binning_df)
    output:
        index=temp("output/binning/minimap2/indexed_contigs/{to_sample}.mmi")
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
        aln=temp("output/binning/minimap2/mapped_reads/{from_sample}.MappedTo.{to_sample}.bam")
    params:
        x=config['params']['minimap2']['x'],
        k=config['params']['minimap2']['k']
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['minimap2_map_reads']
    benchmark:
        "output/benchmarks/map_reads_minimap2/{from_sample}.MappedTo.{to_sample}.benchmark.txt"
    log:
        "output/logs/map_reads_minimap2/{from_sample}.MappedTo.{to_sample}.bowtie.log"
    shell:
        """
        # Map reads against contigs
        minimap2 -a {input.db} {input.reads} -x {params.x} -K {params.k} -t {threads} \
        2> {log} | samtools view -bS - > {output.aln}

        """

rule sort_bam_minimap2:
    """
    Sorts a bam file.
    """
    input:
        aln=rules.map_reads_minimap2.output.aln
    output:
        bam="output/binning/minimap2/mapped_reads/{from_sample}.MappedTo.{to_sample}.sorted.bam"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/sort_bam/minimap2/{from_sample}.MappedTo.{to_sample}.sorted.txt"
    log:
        "output/logs/sort_bam/minimap2/{from_sample}.MappedTo.{to_sample}.sorted.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}

        """
