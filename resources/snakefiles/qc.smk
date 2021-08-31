rule fastqc_pre_trim:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.unit,
                                   wildcards.read)
    output:
        html="output/qc/fastqc/pre_trim/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc/pre_trim/{sample}.{unit}.{read}_fastqc.zip" #  the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    threads:
        config['threads']['fastqc']
    wrapper:
        "0.72.0/bio/fastqc"

rule cutadapt_pe:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.unit,
                                   'R1'),
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.unit,
                                   'R2')
    output:
        fastq1=temp("output/trimmed/{sample}.{unit}.R1.fastq.gz"),
        fastq2=temp("output/trimmed/{sample}.{unit}.R2.fastq.gz"),
        qc="output/logs/cutadapt/{sample}.{unit}.qc.txt"
    params:
        "-a {} {}".format(config["params"]["cutadapt"]['adapter'],
                          config["params"]["cutadapt"]['other'])
    threads:
        config['threads']['cutadapt_pe']
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule fastqc_post_trim:
    input:
        "output/trimmed/{sample}.{unit}.{read}.fastq.gz"
    output:
        html="output/qc/fastqc/post_trim/{sample}.{unit}.{read}.trimmed.html",
        zip="output/qc/fastqc/post_trim/{sample}.{unit}.{read}.trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    threads:
        config['threads']['fastqc_post_trim']
    wrapper:
        "0.72.0/bio/fastqc"


rule merge_units:
    input:
        lambda wildcards: expand("output/trimmed/{sample}.{sequnit}.{read}.fastq.gz",
                                 sample=wildcards.sample,
                                 sequnit=list(units_table.loc[wildcards.sample].index),
                                 read=wildcards.read)
    output:
        temp("output/trimmed/{sample}.combined.{read}.fastq.gz")
    params: ""
    threads: 1
    shell: "cat {input} > {output}"

rule download_NCBI_assembly:
    """
    Downloads a genome by its GenBank assembly accession number.
    """
    output:
        temp(join(config['host_filter']['db_dir'],
                  '{accn}.fa'))
    params:
        accn=config['host_filter']['accn']
    threads: 1
    log: "output/logs/setup/download_NCBI_assembly.{accn}.log"
    conda:
        "../env/bowtie2.yaml"
    shell: "esearch -db assembly -query {params.accn} | \
            elink -target nucleotide -name assembly_nuccore_insdc | \
            efetch -format fasta > {output}"

rule host_bowtie2_build:
    input:
        reference=expand(rules.download_NCBI_assembly.output,
                         accn=config['host_filter']['accn'])
    output:
        multiext(join(config['host_filter']['db_dir'],
                      config['host_filter']['accn']),
                 ".1.bt2",
                 ".2.bt2",
                 ".3.bt2",
                 ".4.bt2",
                 ".rev.1.bt2",
                 ".rev.2.bt2")
    log:
        "output/logs/host_bowtie2_build/bowtie2-build.log"
    conda:
        "../env/bowtie2.yaml"
    params:
        extra="",  # optional parameters,
        indexbase=join(config['host_filter']['db_dir'],
                       config['host_filter']['accn'])
    threads:
        config['threads']['host_filter']
    shell:
        """
        bowtie2-build --threads {threads} {params.extra} \
        {input.reference} {params.indexbase} 2> {log} 1>&2
        """

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
        fastq1="output/trimmed/{sample}.combined.R1.fastq.gz",
        fastq2="output/trimmed/{sample}.combined.R2.fastq.gz",
        db=rules.host_bowtie2_build.output
    output:
        nonhost_R1="output/filtered/nonhost/{sample}.1.fastq.gz",
        nonhost_R2="output/filtered/nonhost/{sample}.2.fastq.gz",
        host="output/filtered/host/{sample}.bam",
    params:
        ref=join(config['host_filter']['db_dir'],
                 config['host_filter']['accn'])
    conda:
        "../env/bowtie2.yaml"
    threads:
        config['threads']['host_filter']
    benchmark:
        "output/benchmarks/host_filter/{sample}_benchmark.txt"
    log:
        "output/logs/host_filter/{sample}.bowtie.log"
    shell:
        """
        # Map reads against reference genome
        bowtie2 -p {threads} -x {params.ref} \
          -1 {input.fastq1} -2 {input.fastq2} \
          --un-conc-gz {wildcards.sample}_nonhost \
          2> {log} | samtools view -bS - > {output.host}

        # rename nonhost samples
        mv {wildcards.sample}_nonhost.1 output/filtered/nonhost/{wildcards.sample}.1.fastq.gz
        mv {wildcards.sample}_nonhost.2 output/filtered/nonhost/{wildcards.sample}.2.fastq.gz
        """

rule multiqc:
    input:
        expand("output/qc/fastqc/pre_trim/{units.Index[0]}.{units.Index[1]}.{read}.html",
               units=units_table.itertuples(), read=reads),
        expand("output/logs/cutadapt/{units.Index[0]}.{units.Index[1]}.qc.txt",
               units=units_table.itertuples()),
        expand("output/qc/fastqc/post_trim/{units.Index[0]}.{units.Index[1]}.{read}.trimmed.html",
               units=units_table.itertuples(), read=reads),
        lambda wildcards: expand(rules.host_filter.log,
                                 sample=samples,
                                 read=reads)
    output:
        "output/qc/multiqc/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/multiqc/multiqc.log"
    wrapper:
        "0.72.0/bio/multiqc"
