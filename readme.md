# DMRworkflow: automating differential methylation analysis of DNA methylation sequencing data

This workflow is being developed as part of a [graduate research project](supplementary/project_background.md) to analyze whole genome bisulfite sequencing (WGBS) data from the Sequence Read Archive (SRA). 
DMRworkflow uses Snakemake workflow management in conjunction with Slurm workload manager to improve the automation and reproducibility of genomic data analysis on an HPC.

### Analysis steps included in pipeline

The snakemake workflow includes a subworkflow component to download fastq files from SRA, run QC checks, trimming, alignment, and methylation calling. 
For human WGBS data analysis, these subworkflow steps require compute resources that are generally not available locally so the subworkflow has been configured to run on an HPC cluster.
The output files from the subworkflow are portable `.beta` files that contain the methylation beta values for each CpG site in the genome.
The primary snakemake workflow requires far less computing resources and can be carried out locally to identify differentially methylated regions (DMRs) associated with a group of samples relative to the background/control group.
Differentially methylated genes (DMGs) are identified as the gene associated with promoters within DMRs.

#### Simplified overview of pipeline

<img src="supplementary/imgs/dmr_workflow_simplified.png" width=75% height=75%>

### Setup

*This pipeline is a work and progress and instructions for setup will be updated as the pipeline is developed.*

##### Dependencies
* python 3.11.4
* conda 23.5.0
    * Automatically manages all Snakemake rule dependencies as specified in `.yaml` files.
    * The `mamba` package manager could be used in place of `conda` for faster dependency resolution.
* snakemake version 7.28.3
    * Create a snakemake environment with `conda create -n snakemake -c conda-forge -c bioconda snakemake`
* wgbstools version 0.2.0
    * This package must be set up, compiled, and in path following the instructions in the [wgbstools repository](https://github.com/nloyfer/wgbs_tools).
    * You do not need to initialize a reference genome with wgbstools, as the pipeline will do this for you.

*More set up and usage instructions to come...*