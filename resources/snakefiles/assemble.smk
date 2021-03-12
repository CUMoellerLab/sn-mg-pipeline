rule metaspades_assembly:
    """

    Performs a metagenomic assembly on a sample using MetaSPAdes.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        contigs="output/metaspades/{sample}.contigs.fasta",
        temp_dir=temp(directory("output/{sample}_temp"))
    conda:
        "resources/envs/spades.yaml"
    threads:
        config['threads']['assembly']
    benchmark:
        "output/benchmarks/metaspades/{sample}_benchmark.txt"
    log:
        "output/logs/metaspades/{sample}.metaspades.log",
    shell:
        """
        # Make temporary output directory
        mkdir -p {output.temp_dir}

        # run the metaspades assmebly
        metaspades.py --threads {threads} \
          -o {output.temp_dir}/ \
          --pe1-1 {input.fastq1} \
          --pe1-2 {input.fastq2} \
          2> {log} 1>&2

        # move and rename the contigs file into a permanent directory
        mv {output.temp_dir}/contigs.fasta {output.contigs}

        """

        