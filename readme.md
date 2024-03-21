# DMRworkflow: automating differential methylation analysis of DNA methylation sequencing data

This workflow is being developed as part of a graduate research project to analyze whole genome bisulfite sequencing (WGBS) data from the Sequence Read Archive (SRA). 
DMRworkflow uses Snakemake workflow management in conjunction with Slurm workload manager to improve the automation and reproducibility of genomic data analysis on an HPC. 

### Pipeline

The snakemake workflow includes a subworkflow component to download fastq files from SRA, run QC checks, trimming, alignment, and methylation calling. 
For human WGBS data analysis, these subworkflow steps require compute resources that are generally not available locally so the subworkflow has been configured to run on an HPC cluster.
The output files from the subworkflow are portable `.beta` files that contain the methylation beta values for each CpG site in the genome. 
The primary snakemake workflow requires far less computing resources and can be carried out locally to identify differentially methylated regions (DMRs) associated with a group of samples relative to the background/control group.
Differentially methylated genes (DMGs) are identified as the gene associated with promoters within DMRs.

### Dependencies:
* Python
* Snakemake
* [wgbstools](https://github.com/nloyfer/wgbs_tools): must be set up, compiled, and in path
* Conda: used by snakemake workflow to manage dependencies for environments as specified in `.yaml` files

## Identifying differentially methylated regions (DMRs) associated with colorectal cancer

The thesis project associated with this workflow seeks to aggregate genome wide methylation data for colorectal cancer patients and healthy patients to identify common differentially methylated regions (DMRs) between sample groups.


### :triangular_flag_on_post: Project Background

#### Methylation patterns observed in cancers:

Epigenome-wide association studies (EWAS) investigate relationships between epigenetic modifications across the entire genome and a particular condition. 

Differentially methylated regions (DMRs) associated with a condition are identified by comparing methylation in two groups.


![Screen Shot 2023-12-12 at 12.10.34 PM](https://hackmd.io/_uploads/r19wRN8UT.png)


DNAm patterns observed in cancer cells:

| Genomic feature    | Change           | Impact of change           |
| ------------------ | ---------------- | -------------------------- |
| Intergenic repeats | Hypomethylation  | Genomic instability        |
| Gene promoters     | Hypomethylation  | Gene reactivation          |
| Gene promoters     | Hypermethylation | Gene silencing             |
| Enhancer           | Hypermethylation | Reduces gene transcription |
| CTCF binding sites | Both             | Genomic instability        |


Methylation arrays, despite covering ~2% of potential methylated regions, are commonly used due to cost considerations. 

Whole genome methylation sequencing, is limited in usage due to high costs, often restricting analysis to one cancer tissue sample per study. 

There is value in aggregating WGBS data from various studies and analysing as a collective.