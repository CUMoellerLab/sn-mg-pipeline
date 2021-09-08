import pandas as pd
from os.path import join

configfile: "config.yaml"

samples_fp = config['samples']

units_fp = config['units']

reads = config['reads']


sample_table = pd.read_csv(samples_fp, sep='\t', header=0)
sample_table.set_index('Sample', inplace=True)

units_table = pd.read_csv(units_fp, sep='\t', header=0)
units_table.set_index(['Sample', 'Unit'], inplace=True)

samples = sample_table.index
units = units_table.index

def get_read(sample, unit, read):
    return(units_table.loc[(sample, unit), read])

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/prototype_selection.smk"
include: "resources/snakefiles/profile.smk"

rule all:
    input:
        "output/qc/multiqc/multiqc.html",
        "output/assemble/multiqc_assemble/multiqc.html",
        "output/prototype_selection/sourmash_plot",
        "output/prototype_selection/prototype_selection/selected_prototypes.yaml",
        "output/profile/metaphlan/merged_abundance_table.txt"
