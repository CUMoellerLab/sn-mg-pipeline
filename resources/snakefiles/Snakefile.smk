import pandas as pd

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

rule all:
    input:
        expand("output/qc/fastqc/{units.Index[0]}.{units.Index[1]}.{read}.html",
        	   units=units_table.itertuples(), read=reads)
