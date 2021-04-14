rule metaspades_assembly:
    """

    Performs a metagenomic assembly on a sample using MetaSPAdes.

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

rule megahit_assembly:
    """

    Performs a metagenomic assembly on a sample using MetaSPAdes.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        contigs="output/assemble/megahit/{sample}.contigs.fasta",
        temp_dir=temp(directory("output/{sample}_temp"))
    conda:
        "../env/assemble.yaml"
    threads:
        config['threads']['megahit']
    benchmark:
        "output/benchmarks/megahit/{sample}_benchmark.txt"
    log:
        "output/logs/megahit/{sample}.megahit.log"
    resources:
        mem_mb=config['mem_mb']['megahit']
    shell:
        """
        # run the metaspades assembly
        megahit -t {threads} \
          -o {output.temp_dir}/ \
          --memory $(({resources.mem_mb}*1024*1024)) \
          -1 {input.fastq1} \
          -2 {input.fastq2} \
          2> {log} 1>&2

        # move and rename the contigs file into a permanent directory
        mv {output.temp_dir}/final.contigs.fa {output.contigs}

        """

rule quast:
    """
    Does an evaluation of assembly quality with Quast
    """
    input:
        lambda wildcards: expand("output/assemble/{assembler}/{sample}.contigs.fasta",
                                 assembler=config['assemblers'],
                                 sample=wildcards.sample)
    output:
        report="output/assemble/quast/{sample}/report.txt",
        outdir=directory("output/assemble/quast/{sample}")
    threads:
        1
    log:
        "output/logs/quast/{sample}.quast.log"
    conda:
        "../env/assemble.yaml"
    benchmark:
        "output/benchmarks/quast/{sample}_benchmark.txt"
    shell:
        """
        quast.py \
          -o {output.outdir} \
          -t {threads} \
          {input}
        """


rule multiqc_assemble:
    input:
        lambda wildcards: expand("output/assemble/quast/{sample}/report.txt",
                                 sample=samples)
    output:
        "output/assemble/multiqc/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/multiqc/multiqc_assemble.log"
    wrapper:
        "0.72.0/bio/multiqc"
