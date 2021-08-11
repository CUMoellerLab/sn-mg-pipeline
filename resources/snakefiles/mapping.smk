from os.path import splitext

rule index_contigs:
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

rule map_reads:
    """
    Performs host read filtering on paired end data using Bowtie and Samtools/
    BEDtools.

    Also requires an indexed reference (path specified in config).

    First, uses Bowtie output piped through Samtools to only retain read pairs
    that are never mapped (either concordantly or just singly) to the indexed
    reference genome. Fastqs from this are gzipped into matched forward and
    reverse pairs.

    Unpaired forward and reverse reads are simply run through Bowtie and
    non-mapping gzipped reads output.

    All piped output first written to localscratch to avoid tying up filesystem.
    """
    input:
        reads = lambda wildcards: expand("output/filtered/nonhost/{sample}.{read}.fastq.gz",
                                         sample=wildcards.from_sample,
                                         read=['1','2']),
        db=rules.index_contigs.output
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
        "output/benchmarks/map_reads/{from_sample}.{to_sample}.benchmark.txt"
    log:
        "output/logs/map_reads/{from_sample}.{to_sample}.bowtie.log"
    shell:
        """
        # Map reads against reference genome
        {params.bt2_command} {params.extra} -p {threads} -x {params.ref} \
          -1 {input.reads[0]} -2 {input.reads[1]} \
          2> {log} | samtools view -bS - > {output.aln}

        """

rule sort_bam:
    """
    Sorts a bam file.
    """
    input:
        aln="output/binning/mapped_reads/{from_sample}.{to_sample}.bam"
    output:
        bam="output/binning/mapped_reads/{from_sample}.{to_sample}.sorted.bam"
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['sort_bam']
    benchmark:
        "output/benchmarks/sort_bam/{from_sample}.{to_sample}.sorted.txt"
    log:
        "output/logs/sort_bam/{from_sample}.{to_sample}.sort.log"
    shell:
        """
        samtools sort -o {output.bam} -@ {threads} {input.aln} 2> {log}

        """
