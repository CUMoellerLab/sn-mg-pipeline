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

binning_fp = config['binning']

binning_df = pd.read_csv(binning_fp, 
                         header=0,
                         index_col=0,
                         sep='\t',
                         na_filter=False)

def get_read(sample, unit, read):
    return(units_table.loc[(sample, unit), read])

def parse_groups(group_series):
    groups = {}
    for sample, grps in group_series.iteritems():
        if not grps:
            continue
        grp_list = grps.split(',')
        for grp in grp_list:
            if grp not in groups:
                groups[grp] = [sample]
            else:
                groups[grp].append(sample)
    return(groups)

def make_pairings(from_grp, to_grp):
    if from_grp.keys() != to_grp.keys():
        raise ValueError('Not all keys in both from and to groups!')
    
    to_samples = []
    from_samples = []
    for grp in from_grp.keys():
        f = from_grp[grp]
        t = to_grp[grp]
        
        for i in f:
            for  j in t:
                from_samples.append(i)
                to_samples.append(j)
    
    return(from_samples, to_samples)

print(binning_df)

to_groups = parse_groups(binning_df['To_Groups'])
from_groups = parse_groups(binning_df['From_Groups'])
from_samples, to_samples = make_pairings(from_groups, to_groups)

print(to_groups)
print(from_groups)
print(from_samples)
print(to_samples)



def get_contigs(sample, binning_df):
    return(binning_df.loc[sample, 'Contigs'])


include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/mapping.smk"

rule map_all:
    input:
        expand("output/binning/mapped_reads/{from_sample}.{to_sample}.sortd.bam",
               from_sample=from_samples,
               to_sample=to_samples)

# rule map_pair:
#     input: 
#         contigs = lambda wildcards: get_contigs(wildcards.to_sample, binning_df),
#         reads = lambda wildcards: expand(rules.merge_units.output,
#                                          sample=wildcards.from_sample,
#                                          read=['R1','R2'])
#     output:
#         "output/binning/mapped_reads/{from_sample}.{to_sample}.bam"
#     run:
#         contig_fp = get_contigs(wildcards.to_sample, binning_df)
#         print("binning_df: \n%s" % binning_df)
#         print("wildcard: %s" % wildcards.to_sample)
#         print("contigs: %s" % contig_fp)
