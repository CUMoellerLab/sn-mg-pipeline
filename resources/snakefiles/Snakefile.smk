import pandas as pd

configfile: "config.yaml"

samples_fp = config['samples']

sample_table = pd.read_csv(samples_fp, sep='\t', header=0)
sample_table.set_index('Sample_ID', inplace=True)

samples=sample_table.index
reads=['R1', 'R2']

def get_read(sample, read):
    return(sample_table.loc[sample, read])

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"

rule all:
    input:
        "output/qc/multiqc/multiqc.html",
        lambda wildcards: expand("output/metaspades/{sample}.contigs.fasta",
                                 sample=samples)
