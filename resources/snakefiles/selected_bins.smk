
def get_bam_list(sample, mapper, contig_pairings):
    fps = expand("output/mapping/{mapper}/sorted_bams/{contig_pairings}_Mapped_To_{sample}.bam",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fps)

def get_index_list(sample, mapper, contig_pairings):
    fps = expand("output/mapping/{mapper}/sorted_bams/{contig_pairings}_Mapped_To_{sample}.bam.bai",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fps)

rule metabat2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = output/binning/metabat2/minimap2/run_metabat2/
        bams = lambda wildcards: get_bam_list(wildcards.contig_sample, wildcards.mapper, contig_pairings)
        bins = lambda wildcards: expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/{contig_sample}_bin.100.fa",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        coverage_table="output/binning/metabat2/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/metabat2/{mapper}/make_metabat2_coverage_table/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/metabat2/{mapper}/make_metabat2_coverage_table/{contig_sample}.log"
    shell:
        """
            Fasta_to_Scaffolds2Bin.sh \
            -i binning/metabat/ \
            -e fa > binning/das_tool/metabat_scaffolds2bin.tsv

            jgi_summarize_bam_contig_depths --outputDepth {output.coverage_table} {input.bams} 2> {log} 1>&2
        """
