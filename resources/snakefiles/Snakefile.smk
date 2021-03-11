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
        "output/metaspades/S22207_S103.contigs.fasta",
        "output/metaspades/S22205_S104.contigs.fasta",

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
        host_R1="output/filtered/host/{sample}.1.fastq.gz",
        host_R2="output/filtered/host/{sample}.2.fastq.gz",
        temp_dir=temp(directory("output/filtered/{sample}_temp"))
    params:
        ref=config['host_reference']
    conda:
        "resources/envs/bowtie2.yaml"
    threads:
        config['threads']['host_filter']
    log:
        bowtie = "output/logs/bowtie2/{sample}.bowtie.log",
        mapped = "output/logs/bowtie2/{sample}.mapped.log",
        unmapped = "output/logs/bowtie2/{sample}.unmapped.log"
    shell:
        """
        # Make temporary output directory
        mkdir -p {output.temp_dir}

        # Map reads against reference genome
        bowtie2 -p {threads} -x {params.ref} --very-sensitive \
          -1 {input.fastq1} -2 {input.fastq2} \
          2> {log.bowtie} > {output.temp_dir}/{wildcards.sample}.alignment.sam

        # separate all read pairs
        # that map at least once to the reference, even discordantly.
        cat {output.temp_dir}/{wildcards.sample}.alignment.sam | samtools view -f 12 -F 256 -b \
          -o {output.temp_dir}/{wildcards.sample}.unsorted.unmapped.bam \
          2> {log.unmapped}

        # separate all read pairs that map to host genome
        cat {output.temp_dir}/{wildcards.sample}.alignment.sam | samtools view -f 2 -F 256 -b \
          -o {output.temp_dir}/{wildcards.sample}.unsorted.mapped.bam \
          2> {log.mapped}

        # Sort the resulting unmapped alignment
        samtools sort -T {output.temp_dir}/{wildcards.sample}.unmapped \
          -@ {threads} -n \
          -o {output.temp_dir}/{wildcards.sample}.unmapped.bam \
          {output.temp_dir}/{wildcards.sample}.unsorted.unmapped.bam \
          2>> {log.unmapped}

        # Sort the resulting mapped alignment
        samtools sort -T {output.temp_dir}/{wildcards.sample}.mapped \
          -@ {threads} -n \
          -o {output.temp_dir}/{wildcards.sample}.mapped.bam \
          {output.temp_dir}/{wildcards.sample}.unsorted.mapped.bam \
          2>> {log.mapped}

        # Convert sorted unmapped alignment to fastq format
        bedtools bamtofastq -i {output.temp_dir}/{wildcards.sample}.unmapped.bam \
          -fq {output.temp_dir}/{wildcards.sample}.R1.trimmed.filtered.fastq \
          -fq2 {output.temp_dir}/{wildcards.sample}.R2.trimmed.filtered.fastq \
          2>> {log.unmapped}

        # Convert sorted mapped alignment to fastq format
        bedtools bamtofastq -i {output.temp_dir}/{wildcards.sample}.mapped.bam \
          -fq {output.temp_dir}/{wildcards.sample}.R1.trimmed.host.fastq \
          -fq2 {output.temp_dir}/{wildcards.sample}.R2.trimmed.host.fastq \
          2>> {log.mapped}

        # zip the filtered fastqs
        pigz -p {threads} \
          -c {output.temp_dir}/{wildcards.sample}.R1.trimmed.filtered.fastq > \
          {output.temp_dir}/{wildcards.sample}.R1.trimmed.filtered.fastq.gz
        pigz -p {threads} \
          -c {output.temp_dir}/{wildcards.sample}.R2.trimmed.filtered.fastq > \
          {output.temp_dir}/{wildcards.sample}.R2.trimmed.filtered.fastq.gz

        # zip the host fastqs
        pigz -p {threads} \
          -c {output.temp_dir}/{wildcards.sample}.R1.trimmed.host.fastq > \
          {output.temp_dir}/{wildcards.sample}.R1.trimmed.host.fastq.gz
        pigz -p {threads} \
          -c {output.temp_dir}/{wildcards.sample}.R2.trimmed.host.fastq > \
          {output.temp_dir}/{wildcards.sample}.R2.trimmed.host.fastq.gz

        # copy the filtered fastqs to final location
        cp {output.temp_dir}/{wildcards.sample}.R1.trimmed.filtered.fastq.gz \
          {output.nonhost_R1}
        cp {output.temp_dir}/{wildcards.sample}.R2.trimmed.filtered.fastq.gz \
          {output.nonhost_R2}

        # copy the host fastqs to final location
        cp {output.temp_dir}/{wildcards.sample}.R1.trimmed.host.fastq.gz \
          {output.host_R1}
        cp {output.temp_dir}/{wildcards.sample}.R2.trimmed.host.fastq.gz \
          {output.host_R2}
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
          --pe1-2 {input.fastq2}

        # move and rename the contigs file into a permanent directory
        mv output/{output.temp_dir}/contigs.fasta output/metaspades/{wildcards.sample}.contigs.fasta

        # move the log file
        mv output/{output.temp_dir}/spades.log output/logs/metaspades/{wildcards.sample}_log.txt

        # remove the temp directory
        # rm -rf {output.temp_dir}

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
        lambda wildcards: expand(rules.host_filter.log.bowtie,
                                 sample=samples)
    output:
        "output/qc/multiqc/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/multiqc/multiqc.log"
    wrapper:
        "0.72.0/bio/multiqc"
