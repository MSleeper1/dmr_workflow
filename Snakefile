

#### import modules ####
import pandas as pd
import workflow_functions as wf

#### assign config ####
configfile: "config.yaml"

#### sample info ####
sample_info = wf.get_sample_info_df(config["samples_tsv"])

subworkflow pre_processing_workflow:
     workdir: 
          config["subworkflow"]["dir"] # path to pre-processing pipeline 
     snakefile:
          config["subworkflow"]["snakefile"] # path to pre-processing snakefile

rule all:
    input:
        pre_processing_workflow(expand("{data_dir}/06_wgbstools_betas/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.{suf}", data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"], sample=sample_info.itertuples())),
        pre_processing_workflow(expand("{rep_dir}/prep_multiqc_data", rep_dir=config["reports_dir"])),
        expand("{data_dir}/06_wgbstools/blocks.bed", data_dir=config["data_dir"]),
        expand("{data_dir}/06_wgbstools/blocks.bed.{suf}", data_dir=config["data_dir"], suf=["gz", "gz.tbi"]),
        expand("{data_dir}/06_wgbstools/segments.tsv", data_dir=config["data_dir"])

# ##### include rules #####
include:"rules/06_wgbstools_segmentation.smk"


