# fastqc each unit
rule fastqc:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.unit,
                                   wildcards.read)
    output:
        html="output/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    threads: 1
    wrapper:
        "0.72.0/bio/fastqc" 


# merge units

# cutadapt

# set up host filter

# host filter

# multiqc