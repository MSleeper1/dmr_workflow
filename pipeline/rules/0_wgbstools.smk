
rule wgbstools_convert_bam_to_beta:
    ''' convert bams into pat and beta files'''
    input:
        bams = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"]), # sambamba_merge output
        ref = expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"])  # wgbstools_init_ref output
         
    output:
        expand("{data_dir}/wgbstools/betas/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.{suf}", data_dir=config["data"]["dir"], suf=["pat.gz", "pat.gz.csi", "beta"])

    log:
        "../pre-processing/logs/rule-logs/wgbstools_convert_bam_to_beta/{ref}/wgbstools_convert_bam_to_beta-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../env/wgbstools.yaml"

    params:
        outdir = expand("{data_dir}/wgbstools/betas", data_dir=config["data"]["dir"]),
        tempdir = expand("{data_dir}/wgbstools/temp", data_dir=config["data"]["dir"]),
        ref_name = config["ref"]["wgbstools_ref_name"]

    shell: 
        """
        mkdir -p {params.outdir}
        mkdir -p {params.tempdir}
        wgbstools bam2pat -f --out_dir {params.outdir} --mbias --genome {params.ref_name} --temp_dir {params.tempdir} {input.bams} > {log} 2>&1
        """




# rule segment_betas:
#     input:
#         betas = expand("{data_dir}/wgbstools/betas/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.beta", data_dir=config["data"]["dir"])
    
#     output:
#         blocks = expand("{data_dir}/wgbstools/blocks.bed", data_dir=config["data"]["dir"])
#         index = expand("{data_dir}/wgbstools/blocks.bed.{suf}", data_dir=config["data"]["dir"], suf=["gz", "gz.tbi"])

#     log:
#         "../pre-processing/logs/rule-logs/wgbstools_segment/{ref}/wgbstools_segment-{ref}-{patient_id}-{group}-{srx_id}-{layout}.err"

#     conda:
#         "../env/wgbstools.yaml"

#     wildcard_constraints:
#         ref = "&&".join(sample_info["ref"].tolist())

#     shell: 
#         """
#         wgbstools segment --betas {input} --min_cpg 3 --max_bp 2000 -o {output.blocks}
#         wgbstools index {output.blocks}
#         """   
        
# ##### left off here
# rule create_table:
#     conda:
#         "dmr.yml"
#     input:
#         blocks="blocks.small.bed.gz",
#         beta="Lung_STL002.small.beta"
#     output:
#         "Lung_STL002-meth-segments.csv"
    
#     shell: 
#         """
#         wgbstools beta_to_table {input.blocks} --betas {input.beta} | column -t >> {output}
#         """

