
# test logging
rule bismark_mapping_se:
    input:
        index = expand("{bismark_idx_dir}/", bismark_idx_dir=config["ref"]["bismark_idx_dir"]),
        read_se = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", data_dir=config["data"]["dir"])
   
    output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bismark_bwt2/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2.bam", data_dir=config["data"]["dir"]),
        report = expand("{rep_dir}/bismark_bwt2/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/bismark_mapping_se/{ref}/bismark_mapping_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../env/bismark.yaml"

    params:
        output_dir = expand("{data_dir}/trimmed/trim_galore/aligned/bismark_bwt2/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data"]["dir"]),
        bismark_idx_dir = expand("{ref_dir}/bismark", ref_dir = config["ref"]["dir"]),
        temp_dir = expand("{data_dir}/tmp", data_dir=config["data"]["dir"]),
        report_dir = expand("{rep_dir}//bismark_bwt2/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),
        accession = "{accession}"

    wildcard_constraints:
        layout = "se"

    shell:
        """
        echo "making output directory {params.output_dir}" > {log}
        mkdir -p {params.output_dir} 2>>{log}
        echo "running bismark on {input.read_se}"
        bismark --bowtie2 --temp_dir {params.temp_dir} --output_dir {params.output_dir} {params.bismark_idx_dir} {input.read_se} >> {log} 2>>{log}
        echo "moving report file to {params.report_dir}" >> {log}
        mkdir -p {params.report_dir} 2>>{log}
        echo " moving {params.output_dir}/{params.accession}_trimmed_bismark_bt2_SE_report.txt to {params.report_dir}" >> {log}
        mv -f -v --target-directory={params.report_dir} {params.output_dir}/{params.accession}_trimmed_bismark_bt2_SE_report.txt 2>>{log}
        """

# bismark mapping for paired end reads
rule bismark_mapping_pe:
    input:
        index = expand("{bismark_idx_dir}/", bismark_idx_dir=config["ref"]["bismark_idx_dir"]),
        r1 = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_trimmed.fq", data_dir=config["data"]["dir"]),
        r2 = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_trimmed.fq", data_dir=config["data"]["dir"])
   
    output:
        expand("{data_dir}/trimmed/trim_galore/aligned/bismark_bwt2/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}{suf}", data_dir=config["data"]["dir"], suf=["_trimmed_bismark_bt2.bam","_trimmed_bismark_bt2_SE_report.txt"])
    
    log:
        "../pre-processing/logs/rule-logs/bismark_mapping_se/{ref}/bismark_mapping_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../env/bismark.yaml"

    params:
        output_dir = expand("{data_dir}/trimmed/trim_galore/aligned/bismark_bwt2/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data"]["dir"]),
        bismark_idx_dir = expand("{ref_dir}/bismark", ref_dir = config["ref"]["dir"]),
        temp_dir = expand("{data_dir}/tmp", data_dir=config["data"]["dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        mkdir -p {params.output_dir}
        bismark --bowtie2 --temp_dir {params.temp_dir} --output_dir {params.output_dir} {params.bismark_idx_dir} -1 {input.r1} -2 {input.r2} > {log} 2>&1
        """