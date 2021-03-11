# Moeller Lab Metagenomics Processing Pipeline

## Overview
Snakemake pipeline for basic processing of metagenomic data from the lab.

Module 1: **sn-mg-QC**

Inputs:
  - Directory of raw fastq files from samples
  - `File_Paths.txt` document, with 4 columns: `Sample_ID`, `R1`, `R2`, and `Sequencing_Run`
    - Files will be renamed with `Sample_ID`
    - `R1` and `R2` should indicate the full path to the forward and reverse reads, respectively
    - `Sequencing_Run` string information may be used to merge samples with the same `Sample_ID` identifier.
  - `config.yaml` file that defines all of the parameters of the pipeline
    - be sure to specify the NCBI Assembly ID for the correct host genome. If not present in `resources/db/bt2/`, the bowtie2 index of this host genome will be automatically downloaded.

Outputs:
  - Directory of quality filtered fastq files that have been depleted of host reads (`output/nonhost/`)
  - Directory of .BAM alignment files containing host reads (`output/host/`)
  - Directory of MinHash profiles of each sample (`output/sourmash/`)
  - skbio Distance Matrix (`output/dist_mat.txt`)
  - List of prototypical samples that best represents the full set of samples (`output/prototypes.txt`)


General Steps:
  1. Quality trim using [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html).
  2. Remove host reads using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
  3. Create MinHash Profiles from each set of Assemblies using [sourmash](https://sourmash.readthedocs.io/en/latest/), and generate a distance matrix of all samples.
  4. Using this distance matrix, select a set of samples that best represents the full sample set using [prototypeSelection](https://github.com/biocore/wol/tree/master/code/prototypeSelection).
  5. Assemble contigs from your quality filtered reads using [MetaSPAdes](https://cab.spbu.ru/software/meta-spades/)

## Installation

First install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) using `conda`.

We recommend first installing and using mamba:

```
$ conda install -c conda-forge mamba
```

Then install snakemake using mamba:

```
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
$ conda activate snakemake
```
Next, clone this github repository
```
$ git clone https://github.com/CUMoellerLab/sn-mg-pipeline.git
cd sn-mg-pipeline
```

Update the `config.yaml` and `resources/config/samples.txt` files
