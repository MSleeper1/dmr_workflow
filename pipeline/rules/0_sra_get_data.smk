# sample_info = pd.read_table(config["samples_tsv"], dtype=str).set_index(["accession"], drop=False)
# accessions = sample_info['accession'].tolist()
# accessions_se = sample_info[sample_info['end_type'] == 'single_end']['accession'].tolist()
# accessions_pe = sample_info[sample_info['end_type'] == 'paired_end']['accession'].tolist()

rule sra_get_data_se:
    output:
        expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.fastq", data_dir = config["data"]["dir"])
        # single_fastq = "/home/msleeper/scratch/data/raw_sequence_files/single_end/{accession}.fastq"
    
    log:
        stdout = "../logs/sra_get_data/sra_get_data-{ref}-{patient_id}-{group}-{srx_id}-{accession}.out",
        stderr = "../logs/sra_get_data/sra_get_data-{ref}-{patient_id}-{group}-{srx_id}-{accession}.err"

    conda:
        "../env/sra-download.yaml"

    params:
        temp_dir = expand("{data_dir}/temp", data_dir = config["data"]["dir"]),
        output_dir = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/", data_dir = config["data"]["dir"])
        # output_dir = expand("{data_dir}/raw_sequence_files/single_end", data_dir = config["data"]["dir"])
        # output_dir = "/home/msleeper/scratch/data/raw_sequence_files/single_end/"

    shell:
        """
        mkdir -p {params.output_dir} > {log.stdout} 2> {log.stderr}
        mkdir -p {params.temp_dir} >> {log.stdout} 2>> {log.stderr}
        fasterq-dump --temp {params.temp_dir} -O {params.output_dir} {wildcards.accession} >> {log.stdout} 2>> {log.stderr}
        """