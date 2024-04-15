
### trim_galore rules ###
# trim_galore is a wrapper around cutadapt and fastqc that trims adapters and low quality bases from fastq files and generates fastqc reports for the trimmed files
# input: raw fastq files
# output: trimmed fastq files and fastqc reports for the trimmed files

# trim_galore_se: rule to trim single-end fastq files
rule trim_galore_se:
    input:
        expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir=config["data_dir"])

    output:
        trimmed_fq = temporary(expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", data_dir=config["data_dir"])),
        fastqc_reports = expand("{rep_dir}/02_fastqc_post_trim/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], suf=["html", "zip"]),
        trim_reports = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", rep_dir=config["reports_dir"])
        
    log:
        stdout = "../pre-processing/logs/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.out",
        stderr = "../pre-processing/logs/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"

    shadow: 
        "shallow"

    conda:
        "../environment_files/trim_galore.yaml"

    params:
        output_dir=expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data_dir"]),  # trimmed files and reports will be saved in this directory
        temp_fastqc_reports = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", data_dir=config["data_dir"], suf=["html", "zip"]),  # temporary fastqc reports will be saved in this directory by trim galore but moved to a report directory after
        fastqc_rep_dir = expand("{rep_dir}/02_fastqc_post_trim/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),  # fastqc reports will be moved to this directory
        temp_trim_reports = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", data_dir=config["data_dir"]),  # temporary trim reports will be saved in this directory by trim galore but moved to a report directory after
        trim_rep_dir = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),  # trim reports will be moved to this directory
        user_args=config["prep_args"]["trim_galore_se"]  # user args can be adjusted in the config file

    wildcard_constraints:
        layout="se"

    shell:
        """
        mkdir -p {params.output_dir} > {log.stdout} 2> {log.stderr}
        echo "trimming with: trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input}" >> {log.stdout} 2>> {log.stderr}
        trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input} >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.fastqc_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.fastqc_rep_dir} {params.temp_fastqc_reports} >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.trim_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.trim_rep_dir} {params.temp_trim_reports} >> {log.stdout} 2>> {log.stderr}
        """

# trim_galore_pe: rule to trim paired-end fastq files
rule trim_galore_pe:
    input:
        r1 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", data_dir=config["data_dir"]),
        r2 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", data_dir=config["data_dir"])

    output:
        trimmed_fq = temporary(expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_trimmed.fq", data_dir=config["data_dir"], read=["1", "2"])),
        fastqc_reports = expand("{rep_dir}/02_fastqc_post_trim/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], read=["1", "2"], suf=["html", "zip"]),
        trim_reports = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], read=["1", "2"])

    log:
        stdout = "/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_pe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.out",
        stderr = "/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_pe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"


    shadow: 
        "shallow"

    conda:
        "../environment_files/trim_galore.yaml"

    params:
        output_dir=expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data_dir"]),  # trimmed files and reports will be saved in this directory
        temp_fastqc_reports = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_trimmed_fastqc.{suf}", data_dir=config["data_dir"], read=["1", "2"], suf=["html", "zip"]),  # temporary fastqc reports will be saved in this directory by trim galore but moved to a report directory after
        fastqc_rep_dir = expand("{rep_dir}/02_fastqc_post_trim/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),  # fastqc reports will be moved to this directory
        temp_trim_reports = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}.fastq_trimming_report.txt", data_dir=config["data_dir"], read=["1", "2"]),  # temporary trim reports will be saved in this directory by trim galore but moved to a report directory after
        trim_rep_dir = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),  # trim reports will be moved to this directory
        user_args=config["prep_args"]["trim_galore_se"]  # user args can be adjusted in the config file

    wildcard_constraints:
        layout="pe"

    shell:
        """
        mkdir -p {params.output_dir} > {log.stdout} 2> {log.stderr}
        trim_galore {params.user_args} --paired --fastqc --output_dir {params.output_dir} {input.r1} {input.r2} >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.fastqc_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.fastqc_rep_dir} {params.temp_fastqc_reports} >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.trim_rep_dir} >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.trim_rep_dir} {params.temp_trim_reports} >> {log.stdout} 2>> {log.stderr}
        """