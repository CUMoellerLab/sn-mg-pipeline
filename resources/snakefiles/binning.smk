
rule make_metabat2_coverage_table:
    """
    Uses jgi_summarize_bam_contig_depths to generate a depth.txt file.
    """
    input:
#        bams = lambda wildcards: expand("output/binning/{mapper}/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.sorted.bam",
#                        mapper=config['mappers'],
#                        read_sample=wildcards.read_sample,
#                        contig_sample=get_contig_list(read_sample, contig_pairings))
#        bams = lambda wildcards: expand(get_bam_list({sample}, config['mappers'], contig_pairings),
#                        sample=wildcards.read_sample)
        bams = get_bam_list({contig_sample}, config['mappers'], contig_pairings)
    output:
        coverage_table="output/binning/metabat2/{contig_sample}_coverage_table.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/metabat2/make_metabat2_depth_file/{contig_sample}_benchmark.txt"
    log:
        "output/logs/metabat2/make_metabat2_depth_file/{contig_sample}.log"
    shell:
        """
            jgi_summarize_bam_contig_depths --outputDepth {output.coverage_table} {input.bams}
        """



rule make_maxbin2_coverage_table:
    """
       Commands to generate a coverage table using `samtools coverage` for input into maxbin2
    """
    input:
        bams = lambda wildcards: expand("output/binning/{mapper}/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.sorted.bam",
                        mapper=config['mappers'],
                        read_sample=read_groups['A'],
                        contig_sample=contig_groups['A'])
    output:
        coverage_table="output/binning/maxbin2/{contig_sample}_coverage_table.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/maxbin2/make_maxbin2_coverage_table/{contig_sample}_benchmark.txt"
    log:
        "output/logs/maxbin2/make_maxbin2_coverage_table/{contig_sample}.log"
    shell:
        """
          samtools coverage {input.bams} | \
          tail -n +2 | \
          sort -k1 | \
          cut -f1,6 > {output.coverage_table}
       """


rule cut_up_fasta:
    """

    Cut up fasta file in non-overlapping or overlapping parts of equal length.
    Optionally creates a BED-file where the cutup contigs are specified in terms
    of the original contigs. This can be used as input to concoct_coverage_table.py.

    """
    input:
        contigs="output/assemble/metaspades/{contig_sample}.contigs.fasta"
    output:
        bed="output/binning/concoct/{contig_sample}_contigs_10K.bed",
        contigs_10K="output/binning/concoct/{contig_sample}_contigs_10K.fa"

    conda:
        "../env/binning.yaml"
    params:
        chunk_size=config['params']['concoct']['chunk_size'],
        overlap_size=config['params']['concoct']['overlap_size']
    benchmark:
        "output/benchmarks/concoct/cut_up_fasta/{contig_sample}_benchmark.txt"
    log:
        "output/logs/concoct/cut_up_fasta/{contig_sample}.log"
    shell:
        """
          python resources/scripts/cut_up_fasta.py {input.contigs} \
          -c {params.chunk_size} \
          -o {params.overlap_size} \
          --merge_last \
          -b {output.bed} > {output.contigs_10K}
        """


# rule make_concoct_coverage_table:
#     """
#
#     Generates table with per sample coverage depth.
#     Assumes the directory "/output/binning/{mapper}/mapped_reads/" contains sorted and indexed bam files where each contig file has has reads mapped against it from the selected prototypes.
#
#     """
#     input:
#         bed=expand("output/binning/concoct/{contig_sample}_contigs_10K.bed",
#                         contig_sample=contig_groups['A']),
#         bams = lambda wildcards: expand("output/binning/{mapper}/mapped_reads/{read_sample}_Mapped_To_{contig_sample}.sorted.bam",
#                         mapper=config['mappers'],
#                         read_sample=read_groups['A'],
#                         contig_sample=contig_groups['A'])
#
#     output:
#         coverage_table="output/binning/concoct/{contig_sample}_coverage_table.txt"
#     conda:
#         "../env/binning.yaml"
#     benchmark:
#         "output/benchmarks/concoct/make_concoct_coverage_table/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/concoct/make_concoct_coverage_table/{contig_sample}.log"
#     shell:
#         """
#           python resources/scripts/concoct_coverage_table.py {input.bed} \
#           {input.bams} > {output.coverage_table}
#         """
