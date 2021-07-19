rule download_metaphlan_db:
    output:
        directory(config['params']['metaphlan']['db_path'])
    conda:
        "../env/metaphlan.yaml"
    log:
        "output/logs/metaphlan/download_metaphlan_db.log"
    benchmark:
        "output/benchmarks/metaphlan/download_metaphlan_db_benchmark.txt"
    shell:
        """
            # Download the database
            metaphlan --install --bowtie2db {output} \
            2> {log} 1>&2

        """


rule run_metaphlan:
    """

    Performs taxonomic profiling using MetaPhlAn3.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2,
        db_path=rules.download_metaphlan_db.output
    output:
        bt2="output/metaphlan/bowtie2s/{sample}.bowtie2.bz2",
        sam="output/metaphlan/sams/{sample}.sam.bz2",
        profile="output/metaphlan/profiles/{sample}_profile.txt",
    conda:
        "../env/metaphlan.yaml"
    threads:
        config['threads']['metaphlan']
    params:
        other=config['params']['metaphlan']['other'],
    benchmark:
        "output/benchmarks/metaphlan/{sample}_benchmark.txt"
    log:
        "output/logs/metaphlan/{sample}.metaphlan.log"
    shell:
        """
        # run metaphlan
        metaphlan {input.fastq1},{input.fastq2} \
        --input_type fastq \
        --nproc {threads} {params.other} \
        --bowtie2db {input.db_path} \
        --bowtie2out {output.bt2}  \
        -s {output.sam}  \
        -o {output.profile} \
        2> {log} 1>&2
        """

rule merge_metaphlan_tables:
    """

    Merges MetaPhlAn3 profiles into a single table

    """
    input:
        expand(rules.run_metaphlan.output.profile,
               sample=samples)
    output:
        merged_abundance_table="output/metaphlan/merged_abundance_table.txt"
    conda:
        "../env/metaphlan.yaml"
    log:
        "output/logs/metaphlan/metaphlan.merged_abundance_table.log"
    shell:
        """

        # merge metaphlan profiles into abundance table
        merge_metaphlan_tables.py {input} \
        -o {output.merged_abundance_table} \
        2> {log} 1>&2

        """
