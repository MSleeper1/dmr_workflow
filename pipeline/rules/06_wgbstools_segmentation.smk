
rule segment_betas:
    input:
        betas = pre_processing_workflow(expand("{data_dir}/06_wgbstools_betas/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.{suf}", data_dir=config["data"]["dir"], suf=["beta"], sample=sample_info.itertuples()))

    output:
        blocks = expand("{data_dir}/06_wgbstools/blocks.bed", data_dir=config["data"]["dir"]),
        index = expand("{data_dir}/06_wgbstools/blocks.bed.{suf}", data_dir=config["data"]["dir"], suf=["gz", "gz.tbi"]),
        table = expand("{data_dir}/06_wgbstools/segments.tsv", data_dir=config["data"]["dir"])

    log:
        "../pre-processing/logs/rule-logs/wgbstools_segment/wgbstools_segment.log"

    conda:
        "../env/wgbstools.yaml"

    shell: 
        """
        wgbstools segment --betas {input.betas} --min_cpg 3 --max_bp 2000 -o {output.blocks}
        wgbstools index {output.blocks}
        wgbstools beta_to_table {output.blocks} --betas {input.betas} | column -t >> {output.table}
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

