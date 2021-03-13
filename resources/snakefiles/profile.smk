rule sourmash_sketch_reads:
    input:
        R1="output/trimmed/{sample}.combined.R1.fastq.gz",
        R2="output/trimmed/{sample}.combined.R2.fastq.gz"
    output:
        "output/sourmash/sketches/{sample}.sig"
    log:
        "output/logs/sourmash/sourmash_sketch_reads.{sample}.log"
    threads: 1
    conda: "../env/sourmash.yaml"
    params:
        k=config['params']['sourmash']['k'],
        scaled=config['params']['sourmash']['scaled'],
        extra=config['params']['sourmash']['extra']
    shell:
        """
        sourmash sketch dna \
        -p k={params.k},scaled={params.scaled}  \
        {params.extra} \
        -o {output} \
        --merge \
        {input} 2> {log} 1>&2
        """

rule sourmash_dm:
    input:
        expand(rules.sourmash_sketch_reads.output,
               sample=samples)
    output:
        "output/sourmash/sourmash.dm"
    log:
        "output/logs/sourmash/sourmash_dm.log"
    threads: 1
    conda: "../env/sourmash.yaml"
    shell:
        """
        sourmash compare \
        --output {output} \
        {input} 2> {log} 1>&2
        """

