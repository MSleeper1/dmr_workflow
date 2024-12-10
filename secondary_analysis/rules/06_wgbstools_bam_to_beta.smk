# note: --mbias is only compatible with pe samples. If I run into issues with mbias plotting error messages I may need to split se and pe samples at this step

rule wgbstools_convert_bam_to_beta_bwa:
    ''' convert bams into pat and beta files'''
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]), # sambamba_merge output
        ref = expand("{root}/{wgbstools_ref_dir}/{fasta}", root = config["root"], wgbstools_ref_dir = config["ref"]["wgbstools_idx_dir"], fasta = config["ref"]["fasta"])
    
    output:
        expand("{root}/{data_dir}/06_wgbstools_betas_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"])
    
    log:
        "logs/secondary_rules/wgbstools_convert_bam_to_beta_bwa/wgbstools_convert_bam_to_beta_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/wgbstools.yaml"

    # shadow:
    #     "shallow"

    params:
        outdir = expand("{root}/{data_dir}/06_wgbstools_betas_bwa/", root = config["root"], data_dir=config["data_dir"]),
        tempdir = expand("{root}/{data_dir}/06_wgbstools_betas_bwa/temp", root = config["root"], data_dir=config["data_dir"]),
        genome_name = config["ref"]["fasta"]

    shell: 
        """
        mkdir -p {params.outdir}
        mkdir -p {params.tempdir}
        echo "Converting bam files {input.bam} to pat and beta files" > {log}
        wgbstools bam2pat -f --out_dir {params.outdir} --mbias --genome {params.genome_name} --temp_dir {params.tempdir} {input.bam} >> {log} 2>&1
        echo "done" >> {log}
        """

rule wgbstools_convert_bam_to_beta_bis:
    ''' convert bams into pat and beta files'''
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]), # sambamba_merge output
        ref = expand("{root}/{wgbstools_ref_dir}/{fasta}", root = config["root"], wgbstools_ref_dir = config["ref"]["wgbstools_idx_dir"], fasta = config["ref"]["fasta"])
    
    output:
        expand("{root}/{data_dir}/06_wgbstools_betas_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.{suf}", root = config["root"], data_dir=config["data_dir"], suf=["pat.gz", "pat.gz.csi", "beta"])
    
    log:
        "logs/secondary_rules/wgbstools_convert_bam_to_beta_bis/wgbstools_convert_bam_to_beta_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/wgbstools.yaml"

    # shadow:
    #     "shallow"

    params:
        outdir = expand("{root}/{data_dir}/06_wgbstools_betas_bis/", root = config["root"], data_dir=config["data_dir"]),
        tempdir = expand("{root}/{data_dir}/06_wgbstools_betas_bis/temp", root = config["root"], data_dir=config["data_dir"]),
        genome_name = config["ref"]["fasta"]

    shell: 
        """
        mkdir -p {params.outdir}
        mkdir -p {params.tempdir}
        echo "Converting bam files {input.bam} to pat and beta files" > {log}
        wgbstools bam2pat -f --out_dir {params.outdir} --mbias --genome {params.genome_name} --temp_dir {params.tempdir} {input.bam} >> {log} 2>&1
        echo "done" >> {log}
        """
