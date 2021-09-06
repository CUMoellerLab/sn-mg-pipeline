rule fastqc_pre_trim:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.unit,
                                   wildcards.read)
    output:
        html="output/qc/fastqc_pre_trim/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc_pre_trim/{sample}.{unit}.{read}_fastqc.zip" #  the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    benchmark:
        "output/benchmarks/qc/fastqc_pre_trim/{sample}.{unit}.{read}_benchmark.txt"
    log:
        "output/logs/qc/fastqc_pre_trim/{sample}.{unit}.{read}.log"
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
        fastq1=temp("output/qc/cutadapt_pe/{sample}.{unit}.R1.fastq.gz"),
        fastq2=temp("output/qc/cutadapt_pe/{sample}.{unit}.R2.fastq.gz"),
        qc="output/logs/qc/cutadapt_pe/{sample}.{unit}.txt"
    params:
        "-a {} {}".format(config["params"]["cutadapt"]['adapter'],
                          config["params"]["cutadapt"]['other'])
    benchmark:
        "output/benchmarks/qc/cutadapt_pe/{sample}.{unit}_benchmark.txt"
    log:
        "output/logs/qc/cutadapt_pe/{sample}.{unit}.log"
    threads:
        config['threads']['cutadapt_pe']
    wrapper:
        "0.17.4/bio/cutadapt/pe"

rule fastqc_post_trim:
    input:
        "output/qc/cutadapt_pe/{sample}.{unit}.{read}.fastq.gz"
    output:
        html="output/qc/fastqc_post_trim/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc_post_trim/{sample}.{unit}.{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    benchmark:
        "output/benchmarks/qc/fastqc_post_trim/{sample}.{unit}.{read}_benchmark.txt"
    log:
        "output/logs/qc/fastqc_post_trim/{sample}.{unit}.{read}.log"
    params: ""
    log:
        "output/logs/qc/fastqc_post_trim/{sample}_{unit}_{read}.log"
    benchmark:
        "output/benchmarks/qc/fastqc_post_trim/{sample}_{unit}_{read}_benchmark.txt"
    threads:
        config['threads']['fastqc_post_trim']
    wrapper:
        "0.72.0/bio/fastqc"

rule merge_units:
    input:
        lambda wildcards: expand("output/qc/cutadapt_pe/{sample}.{sequnit}.{read}.fastq.gz",
                                 sample=wildcards.sample,
                                 sequnit=list(units_table.loc[wildcards.sample].index),
                                 read=wildcards.read)
    output:
        temp("output/qc/merge_units/{sample}.combined.{read}.fastq.gz")
    benchmark:
        "output/benchmarks/qc/merge_units/{sample}.combined.{read}_benchmark.txt"
    log:
        "output/logs/qc/merge_units/{sample}.combined.{read}.log"
    params: ""
    log:
        "output/logs/qc/merge_units/{sample}_combined_{read}.log"
    benchmark:
        "output/benchmarks/qc/merge_units/{sample}_combined_{read}_benchmark.txt"
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
    log:
        "output/logs/qc/download_NCBI_assembly/{accn}.log"
    benchmark:
        "output/benchmarks/qc/download_NCBI_assembly/{accn}_benchmark.txt"
    conda:
        "../env/qc.yaml"
    shell: "esearch -db assembly -query {params.accn} | \
            elink -target nucleotide -name assembly_nuccore_insdc | \
            efetch -format fasta > {output} \
            2>> {log} 1>&2"

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
        "output/logs/host_bowtie2_build/{accn}.log"
    benchmark:
        "output/benchmarks/qc/host_bowtie2_build/{accn}_benchmark.txt"
    conda:
        "../env/qc.yaml"
    params:
        extra="",  # optional parameters
        indexbase=join(config['host_filter']['db_dir'],
                       config['host_filter']['accn'])
    threads:
        config['threads']['host_filter']
    shell:
        """
        bowtie2-build --threads {threads} {params.extra} \
        {input.reference} {params.indexbase}
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
        fastq1="output/qc/merge_units/{sample}.combined.R1.fastq.gz",
        fastq2="output/qc/merge_units/{sample}.combined.R2.fastq.gz",
        db=rules.host_bowtie2_build.output
    output:
        nonhost_R1="output/qc/host_filter/nonhost/{sample}.R1.fastq.gz",
        nonhost_R2="output/qc/host_filter/nonhost/{sample}.R2.fastq.gz",
        host="output/qc/host_filter/host/{sample}.bam",
    params:
        ref=join(config['host_filter']['db_dir'],
                 config['host_filter']['accn'])
    conda:
        "../env/qc.yaml"
    threads:
        config['threads']['host_filter']
    log:
        "output/logs/qc/host_filter/{sample}.log"
    benchmark:
        "output/benchmarks/qc/host_filter/{sample}_benchmark.txt"
    shell:
        """
        # Map reads against reference genome
        bowtie2 -p {threads} -x {params.ref} \
          -1 {input.fastq1} -2 {input.fastq2} \
          --un-conc-gz {wildcards.sample}_nonhost \
          2> {log} | samtools view -bS - > {output.host} 

        # rename nonhost samples
        mv {wildcards.sample}_nonhost.1 output/qc/host_filter/nonhost/{wildcards.sample}.R1.fastq.gz
        mv {wildcards.sample}_nonhost.2 output/qc/host_filter/nonhost/{wildcards.sample}.R2.fastq.gz
        """

rule fastqc_post_host:
    input:
        "output/qc/host_filter/nonhost/{sample}.R1.fastq.gz"
    output:
        html="output/qc/fastqc_post_trim/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc_post_trim/{sample}.{unit}.{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    benchmark:
        "output/benchmarks/qc/fastqc_post_trim/{sample}.{unit}.{read}_benchmark.txt"
    log:
        "output/logs/qc/fastqc_post_trim/{sample}.{unit}.{read}.log"
    params: ""
    log:
        "output/logs/qc/fastqc_post_trim/{sample}_{unit}_{read}.log"
    benchmark:
        "output/benchmarks/qc/fastqc_post_trim/{sample}_{unit}_{read}_benchmark.txt"
    threads:
        config['threads']['fastqc_post_trim']
    wrapper:
        "0.72.0/bio/fastqc"

rule multiqc:
    input:
        expand("output/qc/fastqc_pre_trim/{units.Index[0]}.{units.Index[1]}.{read}.html",
               units=units_table.itertuples(), read=reads),
        expand("output/logs/qc/cutadapt_pe/{units.Index[0]}.{units.Index[1]}.txt",
               units=units_table.itertuples()),
        expand("output/qc/fastqc_post_trim/{units.Index[0]}.{units.Index[1]}.{read}.html",
               units=units_table.itertuples(), read=reads),
        lambda wildcards: expand(rules.host_filter.log,
                                 sample=samples)
    output:
        "output/qc/multiqc/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/qc/multiqc/multiqc.log"
    benchmark:
        "output/benchmarks/qc/multiqc/multiqc_benchmark.txt"
    wrapper:
        "0.72.0/bio/multiqc"
