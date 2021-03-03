# Moeller Lab Metagenomics Processing Pipeline

## Overview
Snakemake pipeline for basic processing of metagenomic data from the lab.

Module 1: **sn-mg-QC**

Inputs:
  - Directory of raw fastq files from samples
  - `File_Paths.txt` document, with 3 columns: `Sample_ID`, `R1`, `R2`
    - Files will be renamed with `Sample_ID`
    - `R1` and `R2` should be the full path to the forward and reverse reads, respectively
  - Host genome in fasta format or bowtie2 index thereof
  - `config.yaml` file that defines all of the parameters of the pipeline

Outputs:
  - Directory of quality filtered fastq files that have been depleted of host reads
  - Directory of quality filtered fastq files containing only host reads
  - Directory of MinHash profiles of each sample

General Steps:
  1. Quality trim using [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)
  2. Remove host reads using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
  3. Assemble genomes in each sample using [metaSPAdes](https://cab.spbu.ru/software/meta-spades/)
  4. Create MinHash Profiles from each set of Assemblies using [sourmash](https://sourmash.readthedocs.io/en/latest/)

## Installation

First install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) using `conda`.

We recommend first installing and using mamba:

```$ conda install -c conda-forge mamba```

Then install snakemake using mamba:

```$ mamba create -c conda-forge -c bioconda -n snakemake snakemake```
```$ conda activate snakemake ```
