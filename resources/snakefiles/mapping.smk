from os.path import splitext

rule index_contigs_bt2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.contig_sample,
                                                binning_df)
    output:
        temp(multiext("output/mapping/bowtie2/indexed_contigs/{contig_sample}",
                      ".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2"))
    log:
        "output/logs/mapping/bowtie2/indexed_contigs/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mapping/bowtie2/indexed_contigs/{contig_sample}_benchmark.txt"
    conda:
        "../env/mapping.yaml"
    params:
        bt2b_command = config['params']['bowtie2']['bt2b_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
        indexbase = "output/mapping/bowtie2/indexed_contigs/{contig_sample}"
    threads:
        config['threads']['bowtie2_build']
    shell:
        """
        {params.bt2b_command} --threads {threads} \
        {input.contigs} {params.indexbase} 2> {log}
        """

rule map_reads_bt2:
    """
    Maps reads to contig files using bowtie2.
    """
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_contigs_bt2.output
    output:
        aln=temp("output/mapping/bowtie2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.bam")
    params:
        ref="output/mapping/bowtie2/mapped_reads/{contig_sample}",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['map_reads']
    benchmark:
        "output/benchmarks/mapping/bowtie2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.benchmark.txt"
    log:
        "output/logs/mapping/bowtie2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
        """
        # Map reads against reference genome
        {params.bt2_command} {params.extra} -p {threads} -x {params.ref} \
          -1 {input.reads[0]} -2 {input.reads[1]} \
          2> {log} | samtools view -bS - > {output.aln}

        """

rule index_contigs_minimap2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.contig_sample,
                                                binning_df)
    output:
        index = temp("output/mapping/minimap2/indexed_contigs/{contig_sample}.mmi")
    log:
        "output/logs/mapping/minimap2/indexed_contigs/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mapping/minimap2/indexed_contigs/{contig_sample}_benchmark.txt"
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['minimap2_index']
    shell:
        """
        minimap2 -d {output.index} {input.contigs} -t {threads} 2> {log}

        """

rule map_reads_minimap2:
    """
    Maps reads to contig files using minimap2.
    """
    input:
        reads = lambda wildcards: expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.read_sample,
                                         read=['R1', 'R2']),
        db=rules.index_contigs_minimap2.output.index
    output:
        aln=temp("output/mapping/minimap2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.bam")
    params:
        x=config['params']['minimap2']['x'],
        k=config['params']['minimap2']['k']
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['minimap2_map_reads']
    benchmark:
        "output/benchmarks/mapping/minimap2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}_benchmark.txt"
    log:
        "output/logs/mapping/minimap2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
        """
        # Map reads against contigs
        minimap2 -a {input.db} {input.reads} -x {params.x} -K {params.k} -t {threads} \
        2> {log} | samtools view -bS - > {output.aln}

        """

rule sort_index_bam:
    """
    Sorts and indexes bam file.
    """
    input:
        aln="output/mapping/{mapper}/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.bam"
    output:
        bam=temp("output/mapping/{mapper}/sorted_bams/{read_sample}_Mapped_To_{contig_sample}.bam"),
        index=temp("output/mapping/{mapper}/sorted_bams/{read_sample}_Mapped_To_{contig_sample}.bam.bai")
    conda:
        "../env/mapping.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/mapping/{mapper}/sort_index_bam/{read_sample}_Mapped_To_{contig_sample}.txt"
    log:
        "output/logs/mapping/{mapper}/sort_index_bam/sort_index_bam/{read_sample}_Mapped_To_{contig_sample}.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}
        samtools index -b -@ {threads} {output.bam} 2>> {log}
        """
