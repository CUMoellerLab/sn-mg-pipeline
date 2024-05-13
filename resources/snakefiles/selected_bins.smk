from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob


rule metabat2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa > {output.scaffolds2bin}
        """


rule maxbin2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fasta > {output.scaffolds2bin}
        """

rule concoct_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa > {output.scaffolds2bin}
        """

# rule concoct_Fasta_to_Scaffolds2Bin:
#     """
#     Uses perl to create a scaffolds2bin.tsv file from a clustering_merged.csv file.
#     """
#     input:
#         merged = lambda wildcards: expand("output/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_clustering_merged.csv",
#                 mapper = config['mappers'],
#                 contig_sample = wildcards.contig_sample)
#     output:
#         scaffolds2bin="output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
#     conda:
#         "../env/selected_bins.yaml"
#     benchmark:
#         "output/benchmarks/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}.log"
#     shell:
#         """
#             perl -pe "s/,/\tconcoct_bins./g;" {input.merged} > {output.scaffolds2bin}
#         """

rule run_DAS_Tool:
    """
    Selects bins using DAS_Tool
    """
    input:
        metabat2 = lambda wildcards: expand("output/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        maxbin2 = lambda wildcards: expand("output/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        concoct = lambda wildcards: expand("output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                    assembler = config['assemblers'],
                    contig_sample = wildcards.contig_sample)
    output:
        out="output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.tsv"
    params:
        basename = "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}",
        search_engine = config['params']['das_tool']['search_engine']
    conda:
        "../env/selected_bins.yaml"
    threads:
        config['threads']['run_DAS_Tool']
    benchmark:
        "output/benchmarks/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}.log"
    shell:
        """
            DAS_Tool \
            --bins {input.metabat2},{input.maxbin2},{input.concoct} \
            --contigs {input.contigs} \
            --outputbasename {params.basename} \
            --labels metabat2,maxbin2,concoct \
            --write_bins \
            --write_bin_evals \
            --threads {threads} \
            --search_engine {params.search_engine}
        """


rule consolidate_DAS_Tool_bins:
    """
    Consolidates and renames bin fastas generated by DAS_Tool into a single folder
    """
    input:
        "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.tsv"
    output:
        done = touch("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done")
    log:
        "output/logs/selected_bins/{mapper}/consolidate_DAS_Tool_bins/{contig_sample}.log"
    run:
        sample = wildcards.contig_sample 
        fasta_dir = join(dirname(input[0]),
                         sample + '_DASTool_bins')
        output_dir = dirname(output.done)

        fasta_files = glob(join(fasta_dir, '*.fa'))

        for file in fasta_files:
            copyfile(file,
                     join(output_dir,
                          sample + '_' + basename(file)))

rule consolidate_DAS_Tool_bins_all:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())


