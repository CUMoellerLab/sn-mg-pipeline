rule kraken2:
    """

    profiles reads using kraken2

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        contigs="output/assemble/metaspades/{sample}.contigs.fasta",
        temp_dir=temp(directory("output/{sample}_temp"))
    conda:
        "../env/assemble.yaml"

    threads:
        config['threads']['spades']
    benchmark:
        "output/benchmarks/metaspades/{sample}_benchmark.txt"
    log:
        "output/logs/metaspades/{sample}.metaspades.log"
    resources:
        mem_mb=config['mem_mb']['spades']
    shell:
        """
        # Make temporary output directory
        mkdir -p {output.temp_dir}

        # run the metaspades assembly
        metaspades.py --threads {threads} \
          -o {output.temp_dir}/ \
          --memory $(({resources.mem_mb}/1024)) \
          --pe1-1 {input.fastq1} \
          --pe1-2 {input.fastq2} \
          2> {log} 1>&2

        # move and rename the contigs file into a permanent directory
        mv {output.temp_dir}/contigs.fasta {output.contigs}

        """


rule taxonomy_kraken:
    """
    Runs Kraken with Bracken to construct taxonomic profiles.
    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        report = "output/profile/kraken2/{sample}.report.txt",
        read_out = "output/profile/kraken2/{sample}.output.txt",
        profile = "output/profile/kraken2/{sample}.profile.txt"
    params:
        db = config['params']['kraken2']['db'],
        levels = config['params']['kraken2']['      '],
        map = config['params']['kraken2']['map']
    threads:
        config['threads']['kraken2']
    log:
        taxonomy_dir + "logs/taxonomy_kraken.sample_{sample}.log"
    benchmark:
        "benchmarks/taxonomy/taxonomy_kraken.sample_{sample}.txt"
    shell:
        """
          # get stem file path
          stem={output.report}
          stem=${{stem%.report.txt}}

          # run Kraken to align reads against reference genomes
          kraken2 {input.fastq1} {input.fastq2} \
            --db {params.db} \
            --paired \
            --gzip-compressed \
            --only-classified-output \
            --threads {threads} \
            --report {output.report} \
            --output {output.read_out} \
            2> {log}

          # run Bracken to re-estimate abundance at given rank
          if [[ ! -z {params.levels} ]]
          then
            IFS=',' read -r -a levels <<< "{params.levels}"
            for level in "${{levels[@]}}"
            do
              bracken \
                -d {params.db} \
                -i {output.report} \
                -t 10 \
                -l $(echo $level | head -c 1 | tr a-z A-Z) \
                -o $stem.redist.$level.txt \
                2>> {log} 1>&2
              rm $stem.report_bracken.txt
            done
          fi
          """

rule kraken:
    input:
        expand("output/profile/kraken2/{sample}.report.txt",
               sample=samples)
