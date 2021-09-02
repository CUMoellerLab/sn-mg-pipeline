# Moeller Lab Metagenomics Processing Pipeline

## Overview
Snakemake pipeline for basic processing of metagenomic data from the lab. It accepts raw fastq files of metagenomic data, quality filters it, removes reads that map to the host genome, then builds assemblies of each sample and generates a [sourmash](https://sourmash.readthedocs.io/en/latest/) profile. The current version also generates a taxonomic profile of each sample using [MetaPhlAn3](https://huttenhower.sph.harvard.edu/metaphlan/). Modules that are currently underdevelopment will handle automated binning procedures, as well as strain-level profiling.

## Quick Start Guide

### Install

First, clone this github repository:

```
$ git clone https://github.com/CUMoellerLab/sn-mg-pipeline.git
cd sn-mg-pipeline
```

We recommend installing and using mamba:

```
$ conda install -c conda-forge mamba
```

Then install the snakemake version for this workflow using mamba:

```
$ mamba env create -n snakemake -f resources/env/snakemake.yaml
$ conda activate snakemake
```



### Update Files

Now you can update three files that are located in the `./resources/config` directory.

The first two files are `samples.txt` and `units.txt`. The `samples.txt` file is your basic metadata file, with each row representing a sample in your dataset and each column containing the corresponding information about that sample. The first column should named "Sample" and should contain the name for each sample. Any addition columns are not used at this step.

The `units.txt` file should have only 4 columns, and each row should correspond to a sample found in the `samples.txt` file. The first column, "Sample", should be all or a subset of the "Sample" column in the `samples.txt` file. The second column, "Unit", should denote which analysis block each sample belongs to. In our case, we use sequencing run/lane, but you may use other information based on your experimental design. The third and fourth columns (named "R1" and "R2", respectively) should include the full file paths to the forward and reverse fastq files for that sample.

The last file to update is the the `config.yaml` file. This is where you can select the parameters for each step in the analysis pipeline. Refer to the documentation for each tool individually for more information. Also, be sure to change the NCBI GenBank Accession number to your host of interest.

NOTE: You can select which metagenomic assembler you want to use under the the "assemblers:" header. The current options are [metaSPAdes](https://cab.spbu.ru/software/meta-spades/) and [MEGAHIT](https://github.com/voutcn/megahit) Simply delete the assembler you don't want to use. Otherwise both will run.

### Run the Pipeline

After you have updated the files described above, you can start the pipeline. First run:
```
conda install -n base -c conda-forge mamba
```

Then begin the run using:
```
snakemake --cores 8 --use-conda
```
The first time you run this, it may take longer to set up your conda environment. Be sure to select the appropriate number of cores for your analysis.
