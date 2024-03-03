#!/usr/bin/env  snakemake --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --output={cluster.output} --error={cluster.error} --job-name={cluster.name}" --jobs 20 --use-conda --rerun-incomplete --printshellcmds

#### snakemake subworkflow pipeline ####
# run on farm cluster with:  snakemake -s prep-Snakefile.smk --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-per-task={cluster.cpus-per-task} --output={cluster.output} --error={cluster.error} --job-name={cluster.name}" --jobs 20 --use-conda --rerun-incomplete --printshellcmds --rerun-triggers mtime --dry-run
# runs as a subworkflow for the dmr workflow pipeline
# author: Meghan M. Sleeper
# input: 
#    - tsv file specifying data accession numbers and associated metadata for wgbs data on SRA
#    - config file specifying reference genome and file paths
# output: 
#    - DNA methylation beta values for each sample
#    - multiqc report for qc at various steps of file processing
# overview: 
#    *includes the DMR analyses steps that are disk use and memory intensive. 
#    *results in portable files summarizing QC and DNA methylation.
#    - retrieves and prepares reference genome files
#    - retrieves fastq files from SRA 
#    - trims, aligns, deduplicates, and merges runs by SRX_id
#    - produces portable files with methylation beta values for each sample
#    - quality control reporting on various steps are collapsed into a single multiqc report

#### import modules ####
import pandas as pd
import subworkflow_functions as swf

#### assign config ####
configfile: "../config.yaml"

#### sample info ####
sample_info = swf.get_sample_info_df(config["samples_tsv"])
sample_info_se = sample_info[sample_info['layout'] == 'se']
sample_info_pe = sample_info[sample_info['layout'] == 'pe']

#### default rule ####
rule all:
     input:
          expand("{data_dir}/06_wgbstools_betas/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.{suf}", data_dir=config["data"]["dir"], suf=["pat.gz", "pat.gz.csi", "beta"], sample=sample_info.itertuples()), # wgbstools output
          expand("{rep_dir}/prep_multiqc_data/multiqc.html", rep_dir=config["reports_dir"]),
          expand("{rep_dir}/prep_multiqc_data", rep_dir=config["reports_dir"])

##### include rules #####
include: "../rules/00_rsync_get_ref_genome.smk"
include: "../rules/00_bwameth_index_ref_genome.smk"
include: "../rules/00_bismark_index_ref_genome.smk"
include: "../rules/00_wgbstools_init_ref_genome.smk"
include: "../rules/00_fastq_screen_genome_prep.smk"
include: "../rules/01_sra_get_data.smk"
include: "../rules/01_qc_reporting.smk"
include: "../rules/02_trim_galore.smk"
include: "../rules/03_bwameth_mapping.smk"
include: "../rules/03_bismark_mapping.smk"
include: "../rules/04_sambamba_sort_index_dedup.smk"
include: "../rules/04_qc_reporting_post_dedup.smk"
include: "../rules/05_sambamba_merge.smk"
include: "../rules/05_qc_reporting_post_merge.smk"
include: "../rules/06_wgbstools_bam_to_beta.smk"
include: "../rules/06_multiqc_compile_reports.smk"


#### -- Intermediate files produced by prep-Snakefile that could be added to rule all -- ###
          # ### 00 reference genome file prep ###
          # expand("{ref_dir}/{fasta}.fa.gz", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),  # rsync_get_ref_genome output
          # expand("{ref_dir}/{fasta}.fa", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),  # rsync_get_ref_genome unzipped output
          # expand("{bwa_idx_dir}/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t", "c2t.bwt", "c2t.pac", "c2t.ann", "c2t.amb", "c2t.sa"], bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"]), # bwameth_index_ref output
          # expand("{bismark_idx_dir}", bismark_idx_dir=config["ref"]["bismark_idx_dir"]),  # bismark_index_ref output
          # expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"]),  # wgbstools_init_ref output
          # expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/", genomes_dir=config["genomes_dir"]),
          # expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", genomes_dir=config["genomes_dir"]),
          # expand("{ref_dir}/{gtf}.{suf}", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"], suf = ["gtf.gz", "gtf"]), # gtf file for ref genome
          
          # ### 01 sra_get_data raw sequence fastq outputs ###
          # expand("{data_dir}/01_raw_sequence_files/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # sra_get_data se output
          # expand("{data_dir}/01_raw_sequence_files/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf={"_1.fastq", "_2.fastq"}), # sra_get_data pe R1 and R2 output
          
          # ### 01 QC reports ###
          # # 01 fastqc reports
          # expand("{rep_dir}/01_fastqc/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
          # expand("{rep_dir}/01_fastqc/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
          # # 01 fastq_screen reports
          # expand("{rep_dir}/01_fastq_screen/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
          # expand("{rep_dir}/01_fastq_screen/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
     
          # ### 02 trim_galore output ###
          # # 02 single end sample output: trimmed fastq, trim and fastq reports
          # expand("{data_dir}/02_trimmed_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed.fq", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # trim_galore se fastq output
          # expand("{rep_dir}/02_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
          # expand("{rep_dir}/02_fastqc_post_trim/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
          # # 02 paired end sample output: trimmed fastq, trim report, and fastq reports for R1 and R2
          # expand("{rep_dir}/02_fastqc_post_trim/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
          # expand("{data_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed.fq", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # trim_galore pe reports
          # expand("{rep_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim

          # ### 03 bwameth_mapping output ###
          # # 03 single end sample output: bam and mapping report
          # expand("{data_dir}/03_aligned_bwameth/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed.bam", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # bwameth_mapping se output
          # expand("{rep_dir}/03_bwameth/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bwameth_mapping se report from stderr
          # # 03 paired end sample output: bam and mapping report
          # expand("{data_dir}/03_aligned_bwameth/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed.bam", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples()), # bwameth_mapping pe output
          # expand("{rep_dir}/03_bwameth/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bwameth_mapping pe report from stderr

          # ### 03 bismark_mapping output ###
          # # 03 single end sample output: bam and mapping report
          # expand("{data_dir}/03_aligned_bismark_bwt2/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2.bam", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # mapped bam output
          # expand("{rep_dir}/03_bismark_bwt2/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bismark_mapping se report
          # # 03 paired end sample output: bam and mapping report
          # expand("{data_dir}/03_aligned_bismark_bwt2/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2.bam", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples()), # mapped bam output
          # expand("{rep_dir}/03_bismark_bwt2/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bismark_mapping pe report

          # ### 04 sambamba_sort_index_dedup output ###
          # expand("{data_dir}/04_deduped_sambamba/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.{suf}", data_dir=config["data"]["dir"], sample=sample_info.itertuples(), suf=["bam", "bam.bai"]), # sambamba_dedup output

          # ### 04 QC reports ###
          # # 04 fastqc post-de-dup QC on bam output reports
          # expand("{rep_dir}/04_fastqc_post_dedup/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc se output post-dedup
          # # 04 samtools stats output
          # expand("{rep_dir}/04_samtools_post_dedup/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output

          # ### 05 sambamba merge output ###
          # expand("{data_dir}/05_merged_sambamba/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam", data_dir=config["data"]["dir"], sample=sample_info.itertuples()), # sambamba_merge output

          # ### 05 QC reports ###
          # # 05 fastqc post-merge output reports
          # expand("{rep_dir}/05_fastqc_post_merge/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc post-merge output
          # # 05 featureCounts output
          # expand("{rep_dir}/05_feature_counts/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample=sample_info.itertuples()), # featureCounts output
          # # 05 qualimap output report directories
          # expand("{rep_dir}/05_qualimap/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output
          # # 05 samtools stats output 
          # expand("{rep_dir}/05_samtools_post_merge/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output

