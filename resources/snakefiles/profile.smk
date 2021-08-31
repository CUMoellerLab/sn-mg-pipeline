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
        "output/profile/kraken2/krona/{sample}.report.html"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/kraken2/taxonomy_kraken.sample_{sample}.log"
    shell:
        """
        perl ../scripts/kraken2-translate.pl {input} > {output}
        """

rule kraken:
    input:
        expand("output/profile/kraken2/{sample}.report.txt",
               sample=samples),
        expand("output/profile/krona/{sample}.report.html",
               sample=samples)
