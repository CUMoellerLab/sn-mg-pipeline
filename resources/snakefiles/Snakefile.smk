import pandas as pd

configfile: "config.yaml"

samples_fp = config['samples']

sample_table = pd.read_csv(samples_fp, sep='\t', header=0)
sample_table.set_index('Sample_ID', inplace=True)

samples=sample_table.index
reads=['R1', 'R2']

def get_read(sample, read):
    return(sample_table.loc[sample, read])

rule all:
    input:
        "output/qc/multiqc/multiqc.html",
        ["output/metaspades/{sample}.contigs.fasta".format(sample=sample) for sample in samples]

rule fastqc_pre_trim:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.read)
    output:
        html="output/qc/fastqc/pre_trim/{sample}.{read}.html",
        zip="output/qc/fastqc/pre_trim/{sample}.{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    threads:
        config['threads']['fastqc']
    wrapper:
        "0.72.0/bio/fastqc"

rule cutadapt_pe:
    input:
        lambda wildcards: sample_table.loc[wildcards.sample,
                                      'R1'],
        lambda wildcards: sample_table.loc[wildcards.sample,
                                      'R2']
    output:
        fastq1="output/trimmed/{sample}.R1.fastq.gz",
        fastq2="output/trimmed/{sample}.R2.fastq.gz",
        qc="output/logs/cutadapt/{sample}.qc.txt"
    params:
        "-a {} {}".format(config["params"]["cutadapt"]['adapter'],
                          config["params"]["cutadapt"]['other'])
    threads:
        config['threads']['cutadapt_pe']
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule fastqc_post_trim:
    input:
        "output/trimmed/{sample}.{read}.fastq.gz"
    output:
        html="output/qc/fastqc/post_trim/{sample}.{read}.trimmed.html",
        zip="output/qc/fastqc/post_trim/{sample}.{read}.trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    threads:
        config['threads']['fastqc_post_trim']
    wrapper:
        "0.72.0/bio/fastqc"


rule host_filter:
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
        fastq1=rules.cutadapt_pe.output.fastq1,
        fastq2=rules.cutadapt_pe.output.fastq2
    output:
        nonhost_R1="output/filtered/nonhost/{sample}.1.fastq.gz",
        nonhost_R2="output/filtered/nonhost/{sample}.2.fastq.gz",
        host="output/filtered/host/{sample}.bam",
        temp_dir=temp(directory("output/{sample}_temp"))
    params:
        ref=config['host_reference']
    conda:
        "resources/envs/bowtie2.yaml"
    threads:
        config['threads']['host_filter']
    benchmark:
        "output/benchmarks/bowtie2/{sample}_benchmark.txt"
    log:
        "output/logs/bowtie2/{sample}.bowtie.log"
    shell:
        """
        # Make temporary and permanent output directories
        mkdir -p "output/filtered/host/"
        mkdir -p "output/filtered/nonhost/"
        mkdir -p {output.temp_dir}

        # Map reads against reference genome
        bowtie2 -p {threads} -x {params.ref} \
          -1 {input.fastq1} -2 {input.fastq2} \
          --un-conc-gz {wildcards.sample}_nonhost \
          2> {log} | samtools view -bS - > {output.host}

        # rename nonhost samples
        mv {wildcards.sample}_nonhost.1 output/filtered/nonhost/{wildcards.sample}.1.fastq.gz
        mv {wildcards.sample}_nonhost.2 output/filtered/nonhost/{wildcards.sample}.2.fastq.gz
        """

rule metaspades_assembly:
    """

    Performs a metagenomic assembly on a sample using MetaSPAdes.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        contigs="output/metaspades/{sample}.contigs.fasta",
        temp_dir=temp(directory("output/{sample}_temp"))
    conda:
        "resources/envs/spades.yaml"
    threads:
        config['threads']['assembly']
    benchmark:
        "output/benchmarks/metaspades/{sample}_benchmark.txt"
    log:
        "output/logs/metaspades/{sample}.metaspades.log",
    shell:
        """
        # Make temporary output directory
        mkdir -p {output.temp_dir}

        # Make directory for assemblies and log files
        mkdir -p output/metaspades/
        mkdir -p output/logs/metaspades/

        # run the metaspades assmebly
        metaspades.py --threads {threads} \
          -o output/{output.temp_dir}/ \
          --pe1-1 {input.fastq1} \
          --pe1-2 {input.fastq2} \
          2> {log}

        # move and rename the contigs file into a permanent directory
        mv output/{output.temp_dir}/contigs.fasta output/metaspades/{wildcards.sample}.contigs.fasta

        """

rule multiqc:
    input:
        lambda wildcards: expand(rules.fastqc_pre_trim.output.zip,
                                 sample=samples,
                                 read=reads),
        lambda wildcards: expand(rules.fastqc_post_trim.output.zip,
                                 sample=samples,
                                 read=reads),
        lambda wildcards: expand(rules.cutadapt_pe.output.qc,
                                 sample=samples),
        lambda wildcards: expand(rules.host_filter.log,
                                 sample=samples)
    output:
        "output/qc/multiqc/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/multiqc/multiqc.log"
    wrapper:
        "0.72.0/bio/multiqc"
