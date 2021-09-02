from os.path import splitext

rule index_contigs_bt2:
    input:
        contigs = lambda wildcards: get_contigs(wildcards.contig_sample,
                                                binning_df)
    output:
        temp(multiext("output/binning/bowtie2/indexed_contigs/{contig_sample}",
                      ".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2"))
    log:
        "output/logs/binning/bowtie2-build.{contig_sample}.log"
    conda:
        "../env/bowtie2.yaml"
    params:
        bt2b_command = config['params']['bowtie2']['bt2b_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
        indexbase = "output/binning/bowtie2/indexed_contigs/{contig_sample}"
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
                                         sample=wildcards.read_sample,
                                         read=['1','2']),
        db=rules.index_contigs_bt2.output
    output:
        aln=temp("output/binning/bowtie2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.bam")
    params:
        ref="output/binning/bowtie2/indexed_contigs/{contig_sample}",
        bt2_command = config['params']['bowtie2']['bt2_command'],
        extra = config['params']['bowtie2']['extra'],  # optional parameters
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['map_reads']
    benchmark:
        "output/benchmarks/map_reads_bt2/{read_sample}_Mapped_To_{contig_sample}.benchmark.txt"
    log:
        "output/logs/map_reads_bt2/{read_sample}_Mapped_To_{contig_sample}.bowtie.log"
    shell:
        """
        # Map reads against reference genome
        {params.bt2_command} {params.extra} -p {threads} -x {params.ref} \
          -1 {input.reads[0]} -2 {input.reads[1]} \
          2> {log} | samtools view -bS - > {output.aln}

        """

rule index_contigs_minimap2:
    input:
        # contigs = expand("output/assemble/metaspades/{sample}.contigs.fasta",
        #         sample=samples)
        contigs = "output/assemble/metaspades/{contig_sample}.contigs.fasta"
    output:
        index = temp("output/binning/minimap2/indexed_contigs/{contig_sample}.mmi")
    log:
        "output/logs/binning/minimap2.index.{contig_sample}.log"
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
                                         sample=wildcards.read_sample,
                                         read=['1','2']),
        db=rules.index_contigs_minimap2.output.index
    output:
        aln=temp("output/binning/minimap2/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.bam")
    params:
        x=config['params']['minimap2']['x'],
        k=config['params']['minimap2']['k']
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['minimap2_map_reads']
    benchmark:
        "output/benchmarks/map_reads_minimap2/{read_sample}_Mapped_To_{contig_sample}.benchmark.txt"
    log:
        "output/logs/map_reads_minimap2/{read_sample}_Mapped_To_{contig_sample}.minimap2.log"
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
        aln="output/binning/{mapper}/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.bam"
    output:
        bam="output/binning/{mapper}/sorted_bams/{read_sample}_Mapped_To_{contig_sample}.sorted.bam",
        bai="output/binning/{mapper}/sorted_bams/{read_sample}_Mapped_To_{contig_sample}.sorted.bam.bai"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/sort_bam/{mapper}/{read_sample}_Mapped_To_{contig_sample}.sorted.txt"
    log:
        "output/logs/sort_bam/{mapper}/{read_sample}_Mapped_To_{contig_sample}.sorted.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}
        samtools index -b -@ {threads} {output.bam} 2>> {log}
        """
