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
    for sample, grps in group_series.items():
        if not grps:
            continue
        grp_list = grps.split(',')
        for grp in grp_list:
            if grp not in groups:
                groups[grp] = [sample]
            else:
                groups[grp].append(sample)
    return(groups)

def make_pairings(read_grp, ctg_grp):
    if read_grp.keys() != ctg_grp.keys():
        raise ValueError('Not all keys in both from and to groups!')

    pairings = []
    contig_pairings = {}
    for grp in read_grp.keys():
        r = read_grp[grp]
        c = ctg_grp[grp]

        for i in r:
            for  j in c:
                pairings.append((i, j))
                if j not in contig_pairings:
                    contig_pairings[j] = [i]
                else:
                    contig_pairings[j].append(i)

    return(pairings, contig_pairings)

contig_groups = parse_groups(binning_df['Contig_Groups'])
read_groups = parse_groups(binning_df['Read_Groups'])
pairings, contig_pairings = make_pairings(read_groups, contig_groups)

print('Contig samples: %s' % contig_groups)
print('Read samples: %s' % read_groups)
print('Pairings: %s' % pairings)
print('Contig Pairings: %s' % contig_pairings)

def get_contigs(sample, binning_df):
    return(binning_df.loc[sample, 'Contigs'])

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/mapping.smk"
include: "resources/snakefiles/binning.smk"
include: "resources/snakefiles/selected_bins.smk"


rule select_bins:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())

rule bin_all:
    input:
        expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/",
               mapper=config['mappers'],
               contig_sample=contig_pairings.keys()),
        expand("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/",
               mapper=config['mappers'],
               contig_sample=contig_pairings.keys()),
        expand("output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/",
               mapper=config['mappers'],
               contig_sample=contig_pairings.keys())

rule map_all:
    input:
        expand("output/mapping/{mapper}/sorted_bams/{pairing[0]}_Mapped_To_{pairing[1]}.bam",
               mapper=config['mappers'],
               pairing=pairings)
