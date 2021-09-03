rule metaspades:
    """

    Performs a metagenomic assembly on a sample using MetaSPAdes.

    """
    input:
        fastq1=rules.host_filter.output.nonhost_R1,
        fastq2=rules.host_filter.output.nonhost_R2
    output:
        contigs="output/assemble/metaspades/{sample}.contigs.fasta",
    params:
        temp_dir=directory("output/{sample}_temp/")
    conda:
        "../env/assemble.yaml"
    threads:
        config['threads']['spades']
    benchmark:
        "output/benchmarks/assemble/metaspades/{sample}_benchmark.txt"
    log:
        "output/logs/assemble/metaspades/{sample}.metaspades.log"
    resources:
        mem_mb=config['mem_mb']['spades']
    shell:
        """
        # Make temporary output directory
        mkdir -p {params.temp_dir}

        # run the metaspades assembly
        metaspades.py --threads {threads} \
          -o {params.temp_dir}/ \
          --memory $(({resources.mem_mb}/1024)) \
          --pe1-1 {input.fastq1} \
          --pe1-2 {input.fastq2} \
          2> {log} 1>&2

        # move and rename the contigs file into a permanent directory
        mv {params.temp_dir}/contigs.fasta {output.contigs}
        rm -rf {params.temp_dir}
        """

rule megahit:
    """

    Performs a metagenomic assembly on a sample using MEGAHIT.

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
        "output/benchmarks/assemble/megahit/{sample}_benchmark.txt"
    log:
        "output/logs/assemble/megahit/{sample}.megahit.log"
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
        report="output/assemble/{assembler}/quast/{sample}/report.txt",
    params:
        outdir=directory("output/assemble/{assembler}/quast/{sample}/")
    threads:
        1
    log:
        "output/logs/assemble/{assembler}/quast/{sample}.quast.log"
    conda:
        "../env/assemble.yaml"
    benchmark:
        "output/benchmarks/assemble/{assembler}/quast/{sample}_benchmark.txt"
    shell:
        """
        quast.py \
          -o {params.outdir} \
          -t {threads} \
          {input}
          touch {output.report}
        """

rule multiqc_assemble:
    input:
        lambda wildcards: expand("output/assemble/{assembler}/quast/{sample}/report.txt",
                                 assembler=config['assemblers'],
                                 sample=samples)
    output:
        "output/assemble/{assembler}/multiqc_assemble/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/assemble/{assembler}/multiqc_assemble/multiqc_assemble.log"
    benchmark:
        "output/benchmarks/assemble/{assembler}/multiqc_assemble/multiqc_assemble_benchmark.txt"
    wrapper:
        "0.72.0/bio/multiqc"

rule metaquast:
    """
    Does an evaluation of assembly quality with Quast
    """
    input:
        lambda wildcards: expand("output/assemble/{assembler}/{sample}.contigs.fasta",
                                 assembler=config['assemblers'],
                                 sample=wildcards.sample)
    output:
        report="output/assemble/{assembler}/metaquast/{sample}/report.html",
        outdir=directory("output/assemble/{assembler}/metaquast/{sample}")
    threads:
        config['threads']['metaquast']
    log:
        "output/logs/assemble/{assembler}/metaquast/{sample}_metaquast.log"
    params:
        refs=config['params']['metaquast']['reference_dir']
    conda:
        "../env/assemble.yaml"
    benchmark:
        "output/benchmarks/assemble/{assembler}/metaquast/{sample}_benchmark.txt"
    shell:
        """
        metaquast.py \
          -r {params.refs} \
          -o {output.outdir} \
          -t {threads} \
          {input}
        """

rule multiqc_metaquast:
    input:
        lambda wildcards: expand("output/assemble/{assembler}/metaquast/{sample}/report.html",
                                 assembler=config['assemblers'],
                                 sample=samples)
    output:
        "output/assemble/{assembler}/multiqc_metaquast/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "output/logs/assemble/{assembler}/multiqc_metaquast/multiqc_metaquast.log"
    benchmark:
        "output/benchmarks/assemble/{assembler}/multiqc_metaquast/multiqc_metaquast_benchmark.txt"
    wrapper:
        "0.72.0/bio/multiqc"
