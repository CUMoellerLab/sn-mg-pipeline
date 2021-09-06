rule taxonomy_kraken:
    """
    Runs Kraken with Bracken to construct taxonomic profiles.
    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        report = "output/profile/kraken2/{sample}.report.txt"
    params:
        db = config['params']['kraken2']['db'],
        levels = config['params']['kraken2']['levels'],
        bracken_db = config['params']['kraken2']['bracken-db']
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['kraken2']
    log:
        "output/logs/kraken2/taxonomy_kraken.sample_{sample}.log"
    benchmark:
        "output/benchmarks/kraken2/taxonomy_kraken.sample_{sample}.txt"
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
            --output - \
            2> {log}

          # run Bracken to re-estimate abundance at given rank
          if [[ ! -z {params.levels} ]]
          then
            IFS=',' read -r -a levels <<< "{params.levels}"
            for level in "${{levels[@]}}"
            do
              bracken \
                -d {params.bracken_db} \
                -i {output.report} \
                -t 10 \
                -l $(echo $level | head -c 1 | tr a-z A-Z) \
                -o $stem.redist.$level.txt \
                2>> {log} 1>&2
            done
          fi
          """

rule krona:
    input:
        rules.taxonomy_kraken.output.report
    output:
        "output/profile/krona/{sample}.report.html"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/krona/krona.sample_{sample}.log"
    shell:
        """
        perl resources/scripts/kraken2-translate.pl {input} > {input}.temp
        ktImportText -o {output} {input}.temp
        rm {input}.temp
        """

rule kraken:
    input:
        expand("output/profile/kraken2/{sample}.report.txt",
               sample=samples),
        expand("output/profile/krona/{sample}.report.html",
               sample=samples)

rule download_metaphlan_db:
    output:
        directory(config['params']['metaphlan']['db_path'])
    conda:
        "../env/profile.yaml"
    log:
        "output/logs/metaphlan/download_metaphlan_db.log"
    benchmark:
        "output/benchmarks/metaphlan/download_metaphlan_db_benchmark.txt"
    shell:
        """
        if test -f "{output}/mpa_latest"; then
            touch {output}
            echo "DB already installed at {output}"
        else
            metaphlan --install --bowtie2db {output} \
                2> {log} 1>&2
        fi
        """


rule metaphlan:
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
        "../env/profile.yaml"
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

    Merges MetaPhlAn3 profiles into a single table.

    """
    input:
        expand(rules.run_metaphlan.output.profile,
               sample=samples)
    output:
        merged_abundance_table="output/metaphlan/merged_abundance_table.txt"
    conda:
        "../env/profile.yaml"
    log:
        "output/logs/metaphlan/metaphlan.merged_abundance_table.log"
    shell:
        """
        merge_metaphlan_tables.py {input} \
        -o {output.merged_abundance_table} \
        2> {log} 1>&2
        """
