
rule make_metabat2_coverage_table:
    """
    Uses jgi_summarize_bam_contig_depths to generate a depth.txt file.
    """
    input:
        bams = lambda wildcards: get_bam_list(wildcards.contig_sample, config['mappers'], contig_pairings)
    output:
        coverage_table="output/binning/metabat2/{mapper}/{contig_sample}_coverage_table.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/metabat2/{mapper}/make_metabat2_depth_file/{contig_sample}_benchmark.txt"
    log:
        "output/logs/metabat2/{mapper}/make_metabat2_depth_file/{contig_sample}.log"
    shell:
        """
            jgi_summarize_bam_contig_depths --outputDepth {output.coverage_table} {input.bams} 2> {log}
        """

rule run_metabat2:
    """
    Runs Metabat2:
    MetaBAT2 clusters metagenomic contigs into different "bins", each of which should correspond to a putative genome.

    MetaBAT2 uses nucleotide composition information and source strain abundance (measured by depth-of-coverage by aligning the reads to the contigs) to perform binning.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                assembler = config['assemblers'],
                contig_sample = wildcards.contig_sample),
        coverage_table = lambda wildcards: expand("output/binning/metabat2/{mapper}/{contig_sample}_coverage_table.txt",
                mapper=config['mappers'],
                contig_sample=wildcards.contig_sample)
    output:
        bins = directory("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/")
    params:
        basename = "output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/{contig_sample}_bins",
        extra = config['params']['metabat2']['extra'],  # optional parameters
    threads:
        config['threads']['run_metabat2']
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/metabat2/{mapper}/run_metabat2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/metabat2/{mapper}/run_metabat2/{contig_sample}.log"
    shell:
        """
            metabat2 {params.extra} --numThreads {threads} \
            --inFile {input.contigs} \
            --outFile {params.basename} \
            --abdFile {input.coverage_table} \
            2> {log}
            touch {output.bins}
        """


# rule make_maxbin2_coverage_table:
#     """
#        Commands to generate a coverage table using `samtools coverage` for input into maxbin2
#     """
#     input:
#         bams="output/binning/{mapper}/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.sorted.bam"
#     output:
#         coverage_table="output/binning/maxbin2/{mapper}/{read_sample}_Mapped_To_{contig_sample, [A-Za-z0-9_]+}_coverage.txt"
#     conda:
#         "../env/binning.yaml"
#     benchmark:
#         "output/benchmarks/maxbin2/{mapper}/make_maxbin2_coverage_table/{read_sample}_Mapped_To_{contig_sample, [A-Za-z0-9_]+}_benchmark.txt"
#     log:
#         "output/logs/maxbin2/{mapper}/make_maxbin2_coverage_table/{read_sample}_Mapped_To_{contig_sample, [A-Za-z0-9_]+}.log"
#     shell:
#         """
#           samtools coverage {input.bams} | \
#           tail -n +2 | \
#           sort -k1 | \
#           cut -f1,6 > {output.coverage_table}
#        """
#
# rule make_maxbin2_abund_list:
#     """
#        Combines the file paths from 'make_maxbin2_coverage_table' for MaxBin2
#     """
#     input:
#         lambda wildcards: expand("output/binning/maxbin2/{mapper}/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
#                 mapper = wildcards.mapper,
#                 contig_sample = wildcards.contig_sample,
#                 read_sample = contig_pairings[wildcards.contig_sample])
#     output:
#         abund_list = "output/binning/maxbin2/{mapper}/{contig_sample}_abund_list.txt"
#     benchmark:
#         "output/benchmarks/maxbin2/{mapper}/make_maxbin2_abund_list/{contig_sample}_abund_list_benchmark.txt"
#     log:
#         "output/logs/maxbin2/{mapper}/make_maxbin2_abund_list/{contig_sample}_abund_list.log"
#     run:
#         with open(output.abund_list, 'w') as f:
#             for fp in input:
#                 f.write('%s\n' % fp)
#
#
# rule run_maxbin2:
#     """
#     Runs MaxBin2:
#     MaxBin2 clusters metagenomic contigs (assembled contiguous genome fragments) into different "bins", each of which corresponds to a putative population genome. It uses nucleotide composition information, source strain abundance (measured by depth-of-coverage by aligning the reads to the contigs), and phylogenetic marker genes to perform binning through an Expectation-Maximization (EM) algorithm.
#     """
#     input:
#         contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
#                 assembler = config['assemblers'],
#                 contig_sample = wildcards.contig_sample),
#         abund_list = lambda wildcards: expand("output/binning/maxbin2/{mapper}/{contig_sample}_abund_list.txt",
#                 mapper=config['mappers'],
#                 contig_sample=wildcards.contig_sample)
#     output:
#         bins = "output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}_bins"
#     params:
#         prob = config['params']['maxbin2']['prob_threshold'],  # optional parameters
#         extra = config['params']['maxbin2']['extra']  # optional parameters
#     threads:
#         config['threads']['run_maxbin2']
#     conda:
#         "../env/binning.yaml"
#     benchmark:
#         "output/benchmarks/metabat2/{mapper}/run_maxbin2/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/metabat2/{mapper}/run_maxbin2/{contig_sample}.log"
#     shell:
#         """
#             run_MaxBin.pl -thread {threads} -prob_threshold {params.prob} {params.extra} \
#             -contig {input.contigs} \
#             -abund_list {input.abund_list} \
#             -out {output.bins}
#             2> {log}
#         """
#
# rule cut_up_fasta:
#     """
#     Cut up fasta file in non-overlapping or overlapping parts of equal length.
#     Optionally creates a BED-file where the cutup contigs are specified in terms
#     of the original contigs. This can be used as input to concoct_coverage_table.py.
#     """
#     input:
#         contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
#                 assembler = config['assemblers'],
#                 contig_sample = wildcards.contig_sample)
#     output:
#         bed="output/binning/concoct/{mapper}/{contig_sample}_contigs_10K.bed",
#         contigs_10K="output/binning/concoct/{mapper}/{contig_sample}_contigs_10K.fa"
#     conda:
#         "../env/concoct_linux.yaml"
#     params:
#         chunk_size=config['params']['concoct']['chunk_size'],
#         overlap_size=config['params']['concoct']['overlap_size']
#     benchmark:
#         "output/benchmarks/concoct/{mapper}/cut_up_fasta/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/concoct/{mapper}/cut_up_fasta/{contig_sample}.log"
#     shell:
#         """
#           python resources/scripts/cut_up_fasta.py {input.contigs} \
#           -c {params.chunk_size} \
#           -o {params.overlap_size} \
#           --merge_last \
#           -b {output.bed} > {output.contigs_10K} 2> {log}
#         """
#
#
# rule make_concoct_coverage_table:
#     """
#     Generates table with per sample coverage depth.
#     Assumes the directory "/output/binning/{mapper}/mapped_reads/" contains sorted and indexed bam files where each contig file has has reads mapped against it from the selected prototypes.
#
#     """
#     input:
#         bed="output/binning/concoct/{mapper}/{contig_sample}_contigs_10K.bed",
#         bams = lambda wildcards: get_bam_list(wildcards.contig_sample, config['mappers'], contig_pairings)
#     output:
#         coverage_table="output/binning/concoct/{mapper}/{contig_sample}_coverage_table.txt"
#     conda:
#         "../env/concoct_linux.yaml"
#     benchmark:
#         "output/benchmarks/concoct/{mapper}/make_concoct_coverage_table/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/concoct/{mapper}/make_concoct_coverage_table/{contig_sample}.log"
#     shell:
#         """
#           python resources/scripts/concoct_coverage_table.py {input.bed} \
#           {input.bams} > {output.coverage_table} 2> {log}
#         """
