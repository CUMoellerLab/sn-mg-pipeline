
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
            Fasta_to_Scaffolds2Bin.sh \
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
            Fasta_to_Scaffolds2Bin.sh \
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
            Fasta_to_Scaffolds2Bin.sh \
            -i {input.bins} \
            -e fasta > {output.scaffolds2bin}
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
