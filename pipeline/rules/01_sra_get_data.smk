### SRA download rules ###
# This file contains the rules to download the raw sequence files from the SRA database using the fasterq-dump tool.

# sra_get_data_se: rule to download single-end raw sequence files from the SRA database.
rule sra_get_data_se:
    output:
        temporary(expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir = config["data"]["dir"]))
    
    log:
       "../pre-processing/logs/rule-logs/01_sra_get_data/{ref}/01_sra_get_data-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../env/sra-download.yaml"

    params:
        temp_dir = expand("{data_dir}/temp/01-sra-{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir = config["data"]["dir"]),
        output_dir = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir = config["data"]["dir"])

    wildcard_constraints:
        layout = "se"

    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {params.temp_dir} 
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} > {log} 2>&1
        rm -rf {params.temp_dir}
        """

# sra_get_data_pe: rule to download paired-end raw sequence files from the SRA database.
rule sra_get_data_pe:
    output:
        r1 = temporary(expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", data_dir = config["data"]["dir"])),
        r2 = temporary(expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", data_dir = config["data"]["dir"]))
    
    log:
        "../pre-processing/logs/rule-logs/01_sra_get_data/{ref}/01_sra_get_data-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../env/sra-download.yaml"

    params:
        temp_dir = expand("{data_dir}/temp/01-sra-{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir = config["data"]["dir"]),
        output_dir = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir = config["data"]["dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        mkdir -p {params.output_dir} 
        mkdir -p {params.temp_dir} 
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} > {log} 2>&1
        rm -rf {params.temp_dir}
        """