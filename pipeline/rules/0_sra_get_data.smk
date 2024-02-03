

rule sra_get_data_se:
    output:
        expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.fastq", data_dir = config["data"]["dir"])
    
    log:
       "../logs/sra_get_data/sra_get_data-{ref}-{patient_id}-{group}-{srx_id}-{accession}.log"

    conda:
        "../env/sra-download.yaml"

    params:
        temp_dir = expand("{data_dir}/temp", data_dir = config["data"]["dir"]),
        output_dir = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/", data_dir = config["data"]["dir"])

    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.temp_dir} 
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} > {log} 2>&1
        """

rule sra_get_data_pe:
    output:
        r1 = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_1.fastq", data_dir = config["data"]["dir"]),
        r2 = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_2.fastq", data_dir = config["data"]["dir"])
    
    log:
        "../logs/sra_get_data/sra_get_data-{ref}-{patient_id}-{group}-{srx_id}-{accession}.log"

    conda:
        "../env/sra-download.yaml"

    params:
        temp_dir = expand("{data_dir}/temp", data_dir = config["data"]["dir"]),
        output_dir = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/", data_dir = config["data"]["dir"])

    shell:
        """
        mkdir -p {params.output_dir} 
        mkdir -p {params.temp_dir} 
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} > {log} 2>&1
        """