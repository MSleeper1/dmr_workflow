# note: --mbias is only compatible with pe samples. If I run into issues with mbias plotting error messages I may need to split se and pe samples at this step

rule wgbstools_convert_bam_to_beta:
    ''' convert bams into pat and beta files'''
    input:
        bam = expand("{data_dir}/05_merged_sambamba/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]), # sambamba_merge output
        ref = expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"])  # wgbstools_init_ref output

    output:
        expand("{data_dir}/06_wgbstools_betas/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.{suf}", data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"])
    log:
        "../pre-processing/logs/rule-logs/wgbstools_convert_bam_to_beta/{ref}/wgbstools_convert_bam_to_beta-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../environment_files/wgbstools.yaml"

    params:
        outdir = expand("{data_dir}/06_wgbstools_betas/", data_dir=config["data_dir"]),
        tempdir = expand("{data_dir}/06_wgbstools_betas/temp", data_dir=config["data_dir"]),
        ref_name = config["ref"]["wgbstools_ref_name"]

    shell: 
        """
        mkdir -p {params.outdir}
        mkdir -p {params.tempdir}
        wgbstools bam2pat -f --out_dir {params.outdir} --mbias --genome {params.ref_name} --temp_dir {params.tempdir} {input.bam} > {log} 2>&1
        """

