


rule trim_galore_se:
    input:
        expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.fastq", data_dir=config["data"]["dir"])

    output:
        expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}{suf}", data_dir=config["data"]["dir"], suf=[".fq", ".fastq_trimming_report.txt"]),
        expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_trimmed_fastqc.{suf}", data_dir=config["data"]["dir"], suf=["html", "zip"]) # consider moving to reports directory after debugging

    log:
        "../logs/trim_galore_se/trim_galore_se-{ref}-{patient_id}-{group}-{srx_id}-{accession}.log"

    conda:
        "../env/trim_galore.yaml"

    params:
        output_dir=expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/", data_dir=config["data"]["dir"]),  # trimmed files and reports will be saved in this directory
        user_args=config["prep_args"]["trim_galore_se"]  # user args can be adjusted in the config file

    shell:
        """
        mkdir -p {params.output_dir}
        trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input} 
        """
