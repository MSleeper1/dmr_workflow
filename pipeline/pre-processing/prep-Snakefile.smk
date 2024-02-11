#!/usr/bin/env snakemake --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-{pe.layout}r-task={cluster.cpus-{pe.layout}r-task} --output={cluster.output} --error={cluster.error}" --jobs 11 --use-conda --rerun-incomplete --printshellcmds -s

#### snakemake subworkflow pipeline ####
# run on farm with: snakemake -s prep-Snakefile.smk --cluster-config ../cluster.yaml --cluster "sbatch --parsable --time={cluster.time} --mem={cluster.mem_mb} --nodes={cluster.nodes} --cpus-{pe.layout}r-task={cluster.cpus-{pe.layout}r-task} --output={cluster.output} --error={cluster.error}" --jobs 11 --use-conda --rerun-incomplete --printshellcmds

#### import modules ####
import pandas as pd
import subworkflow_functions as swf

#### assign config ####
configfile: "/home/msleeper/workflows/dmr_workflow/pipeline/config.yaml"

#### sample info ####
sample_info = swf.get_sample_info_df(config["samples_tsv"])
sample_info_se = sample_info[sample_info['layout'] == 'se']
sample_info_pe = sample_info[sample_info['layout'] == 'pe']

# #### make a dictionary with srx_id key : accession list value pairs ####
# srx_acc_dict = swf.make_srx_acc_dict(sample_info)

# #### make dictionary with srx_id key : df rows value pairs for merging info ####
# srx_acc_df_dict = swf.make_srx_df_dict(sample_info)
# srx_merged_dict = swf.merge_rows_by_srx_id(srx_acc_df_dict)
# sample_info_by_srx = swf.make_sample_info_by_srx(sample_info)

#### default rule ####
rule all:
     input:
          ### reference genome files ###
          # expand("{ref_dir}/{fasta}.fa.gz", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),  # rsync_get_ref_genome output
          # expand("{ref_dir}/{fasta}.fa", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),  # rsync_get_ref_genome unzipped output
          # expand("{bwa_idx_dir}/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t", "c2t.bwt", "c2t.pac", "c2t.ann", "c2t.amb", "c2t.sa"], bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"]), # bwameth_index_ref output
          # expand("{bismark_idx_dir}", bismark_idx_dir=config["ref"]["bismark_idx_dir"]),  # bismark_index_ref output
          expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"]),  # wgbstools_init_ref output
          
          ### sra_get_data raw sequence files ###
          expand("{data_dir}/raw_sequence_files/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # sra_get_data se output
          expand("{data_dir}/raw_sequence_files/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf={"_1.fastq", "_2.fastq"}), # sra_get_data pe R1 and R2 output
          
          ### fastqc pre-trim output reports ###
          expand("{rep_dir}/fastqc/pre-trim/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc se output
          expand("{rep_dir}/fastqc/pre-trim/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # fastqc pe R1 and R2 output
          
          ### trim_galore output ###
          expand("{data_dir}/trimmed/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed.fq", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # trim_galore se output
          expand("{data_dir}/trimmed/trim_galore/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed.fq", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # trim_galore pe output
          ### trim_galore output trimming reports ###
          expand("{rep_dir}/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se report
          expand("{rep_dir}/trim_galore/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # moved trim_galore pe reports
          ### trim_galore post-trim fastqc output reports ###
          expand("{rep_dir}/fastqc/post-trim/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # moved fastqc se post-trim output
          expand("{rep_dir}/fastqc/post-trim/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # moved fastqc pe post-trim output
          
          ### bwameth_mapping output ###
          expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed.bam", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # bwameth_mapping se output
          expand("{rep_dir}/bwameth/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bwameth_mapping se report from stderr
          expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed.bam", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples()), # bwameth_mapping pe output
          
          ### bismark_mapping output ###
          expand("{data_dir}/trimmed/trim_galore/aligned/bismark_bwt2/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2.bam", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # bismark_mapping se output
          expand("{rep_dir}/bismark_bwt2/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), 

          ### sambamba_sort_index_dedup output ###
          expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{sample.ref}/{sample.patient_id}/{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.{suf}", data_dir=config["data"]["dir"], sample=sample_info.itertuples(), suf=["bam", "bam.bai"]), # sambamba_dedup output

          ### sambamba flagstat output ###

          ### sambamba merge output ###
          expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam", data_dir=config["data"]["dir"], sample=sample_info.itertuples()) # sambamba_merge output

          ### reports ###

          ### removed reporting in sambamba sort index dedup rule - flagstat will be done separately###
          # expand("{rep_dir}/sambamba/{sample.ref}/{sample.patient_id}/{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}{suf}", rep_dir=config["reports_dir"], sample=sample_info.itertuples(), suf=["_trimmed.bam.flagstat", "_trimmed_sorted_dedup.bam.flagstat"]) # sambamba flagstat output

          # ### sambamba merge outputs ###
          # expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{srx.ref}/{srx.patient_id}/{srx.group}-{srx.srx_id}-{srx.layout}/{srx.group}-{srx.srx_id}_trimmed_bwameth_deduped_merged.bam", data_dir=config["data"]["dir"], srx=sample_info_by_srx.itertuples()) # sambamba_merge output


          ### old input format (no longer needed after testing) ###
          # expand("{ref_dir}/{fasta}.fa.gz", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),  # rsync_get_ref_genome output
          # expand("{ref_dir}/{fasta}.fa", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),  # rsync_get_ref_genome unzipped output
          # expand("{bwa_idx_dir}/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t", "c2t.bwt", "c2t.pac", "c2t.ann", "c2t.amb", "c2t.sa"], bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"]), # bwameth_index_ref output
          # expand("{bismark_idx_dir}", bismark_idx_dir=config["ref"]["bismark_idx_dir"]),  # bismark_index_ref output
          # expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"]),  # wgbstools_init_ref output
          # expand("{data_dir}/raw_sequence_files/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq", data_dir=config["data"]["dir"], se=sample_info_se.itertuples()), # sra_get_data se output
          # expand("{data_dir}/raw_sequence_files/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf={"_1.fastq", "_2.fastq"}), # sra_get_data pe R1 and R2 output
          # expand("{rep_dir}/quality/fastqc/pretrim/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc se output
          # expand("{rep_dir}/quality/fastqc/pretrim/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), suf=["_1_fastqc.html", "_2_fastqc.html", "_1_fastqc.zip", "_2_fastqc.zip"]), # fastqc pe R1 and R2 output
          # expand("{data_dir}/trimmed/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}{suf}", data_dir=config["data"]["dir"], se=sample_info_se.itertuples(), suf=["_trimmed.fq",".fastq_trimming_report.txt"]), # trim_galore se output
          # expand("{data_dir}/trimmed/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}{suf}", data_dir=config["data"]["dir"], se=sample_info_se.itertuples(), suf=["_trimmed_fastqc.html","_trimmed_fastqc.zip"]), # fastqc se post-trim output (run by trim_galore) Consider moving to reports directory after troubleshooting
          # expand("{data_dir}/trimmed/trim_galore/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf=["_1_trimmed.fq", "_2_trimmed.fq",".fastq_trimming_report.txt"]), # trim_galore pe output
          # expand("{data_dir}/trimmed/trim_galore/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf=["_1_trimmed_fastqc.html","_1_trimmed_fastqc.zip", "_2_trimmed_fastqc.html","_2_trimmed_fastqc.zip"]), # fastqc pe post-trim output (run by trim_galore) Consider moving to reports directory after troubleshooting
          # expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}{suf}", data_dir=config["data"]["dir"], se=sample_info_se.itertuples(), suf=["_trimmed.bam"]), # bwameth_mapping se output
          # expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{se.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf=["_trimmed.bam"]), # bwameth_mapping pe output
          # expand("{data_dir}/trimmed/trim_galore/aligned/bismark_bwt2/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}{suf}", data_dir=config["data"]["dir"], se=sample_info_se.itertuples(), suf=["_trimmed_bismark_bt2.bam","_trimmed_bismark_bt2_SE_report.txt"]), # bismark_mapping se output
          # expand("{rep_dir}/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se reports
          # expand("{rep_dir}/fastqc/post-trim/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]) # moved fastqc se post-trim output

          # expand("{rep_dir}/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se reports
          # expand("{rep_dir}/trim_galore/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # moved trim_galore pe reports
          # expand("{rep_dir}/trim_galore/{pe.ref}/{pe.patient_id}/{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}{suf}", data_dir=config["data"]["dir"], pe=sample_info_pe.itertuples(), suf=["_1_trimmed_fastqc.html","_1_trimmed_fastqc.zip", "_2_trimmed_fastqc.html","_2_trimmed_fastqc.zip", ".fastq_trimming_report.txt"]), # fastqc pe post-trim output (run by trim_galore) Consider moving to reports directory after troubleshooting
          # expand("{rep_dir}/trim_galore/{se.ref}/{se.patient_id}/{se.group}-{se.srx_id}-{se.layout}/{se.accession}{suf}", data_dir=config["data"]["dir"], se=sample_info_se.itertuples(), suf=["_trimmed_fastqc.html","_trimmed_fastqc.zip", ".fastq_trimming_report.txt"]), # fastqc se post-trim output (run by trim_galore) Consider moving to reports directory after troubleshooting
          

##### include rules #####
include: "../rules/0_rsync_get_ref_genome.smk"
include: "../rules/0_bwameth_index_ref_genome.smk"
include: "../rules/0_bismark_index_ref_genome.smk"
include: "../rules/0_wgbstools_init_ref_genome.smk"
include: "../rules/0_sra_get_data.smk"
include: "../rules/0_fastqc_pretrim.smk"
include: "../rules/0_trim_galore.smk"
include: "../rules/0_bwameth_mapping.smk"
include: "../rules/0_bismark_mapping.smk"
include: "../rules/0_sambamba_sort_index_dedup.smk"
include: "../rules/0_sambamba_merge.smk"
# include: "../rules/0_reports.smk"


