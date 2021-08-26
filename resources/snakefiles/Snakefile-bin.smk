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

def get_bam_list(sample, mapper, contig_pairings):
    fp = expand("output/binning/{mapper}/mapped_reads/{sample}_Mapped_To_{contig_pairings}.sorted.bam",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fp)

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/mapping.smk"
include: "resources/snakefiles/binning.smk"

rule map_all:
    input:
        expand("output/binning/{mapper}/mapped_reads/{pairing[0]}_Mapped_To_{pairing[1]}.sorted.bam",
                mapper=config['mappers'],
                pairing=pairings),
        # expand("output/binning/metabat2/{mapper}/{contig_sample}_coverage_table.txt",
        #         mapper=config['mappers'],
        #         contig_sample=contig_groups['A']),
        directory(expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample, [A-Za-z0-9_]+}/"),
                mapper=config['mappers'],
                contig_sample=contig_groups['A'])
        # expand("output/binning/maxbin2/{mapper}/{pairing[0]}_Mapped_To_{pairing[1]}_coverage.txt",
        #         mapper=config['mappers'],
        #         pairing=pairings),
        # expand("output/binning/maxbin2/{mapper}/{contig_sample}_abund_list.txt",
        #         mapper=config['mappers'],
        #         contig_sample=contig_groups['A']),
        # expand("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}_bins",
        #         mapper=config['mappers'],
        #         contig_sample=contig_groups['A']),
        # expand("output/binning/concoct/{mapper}/{contig_sample}_coverage_table.txt",
        #         mapper=config['mappers'],
        #         contig_sample=contig_groups['A'])
