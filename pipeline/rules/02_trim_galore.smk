
### trim_galore rules ###
# trim_galore is a wrapper around cutadapt and fastqc that trims adapters and low quality bases from fastq files and generates fastqc reports for the trimmed files
# input: raw fastq files
# output: trimmed fastq files and fastqc reports for the trimmed files

# trim_galore_se: rule to trim single-end fastq files
rule trim_galore_se:
    input:
        expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir=config["data"]["dir"])

    output:
        trimmed_fq = temporary(expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", data_dir=config["data"]["dir"])),
        fastqc_reports = expand("{rep_dir}/02_fastqc_post_trim/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], suf=["html", "zip"]),
        trim_reports = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", rep_dir=config["reports_dir"])
        
    log:
        stdout = "../pre-processing/logs/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.out",
        stderr = "../pre-processing/logs/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"

    conda:
        "../env/trim_galore.yaml"

    params:
        output_dir=expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data"]["dir"]),  # trimmed files and reports will be saved in this directory
        user_args=config["prep_args"]["trim_galore_se"],  # user args can be adjusted in the config file
        fastqc_rep_dir = expand("{rep_dir}/02_post_trim_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),  # fastqc reports will be moved to this directory
        trim_rep_dir = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"])  # trim reports will be moved to this directory

    wildcard_constraints:
        layout="se"

    shell:
        """
        mkdir -p {params.output_dir} > {log.stdout} 2> {log.stderr}
        echo "trimming with: trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input}" >> {log.stdout} 2>> {log.stderr}
        trim_galore {params.user_args} --fastqc --output_dir {params.output_dir} {input} >> {log.stdout} 2>> {log.stderr}
        echo "making directory: {params.fastqc_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.fastqc_rep_dir} >> {log.stdout} 2>> {log.stderr}
        echo "moving fastqc reports to {params.fastqc_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.fastqc_rep_dir} {params.output_dir}/{wildcards.accession}_trimmed_fastqc.html {params.output_dir}/{wildcards.accession}_trimmed_fastqc.zip >> {log.stdout} 2>> {log.stderr}
        echo "making directory: {params.trim_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mkdir -p {params.trim_rep_dir} >> {log.stdout} 2>> {log.stderr}
        echo "moving trim reports to {params.trim_rep_dir}" >> {log.stdout} 2>> {log.stderr}
        mv -f -v --target-directory={params.trim_rep_dir} {params.output_dir}/{wildcards.accession}.fastq_trimming_report.txt >> {log.stdout} 2>> {log.stderr}
        """

# trim_galore_pe: rule to trim paired-end fastq files
rule trim_galore_pe:
    input:
        r1 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", data_dir=config["data"]["dir"]),
        r2 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", data_dir=config["data"]["dir"])

    output:
        trimmed_fq = temporary(expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_trimmed.fq", data_dir=config["data"]["dir"], read=["1", "2"])),
        fastqc_reports = expand("{rep_dir}/02_post_trim_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], read=["1", "2"], suf=["html", "zip"]),
        trim_reports = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], read=["1", "2"])

    log:
        "../pre-processing/logs/rule-logs/02_trim_galore_se/{ref}/02_trim_galore_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../env/trim_galore.yaml"

    params:
        output_dir=expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data"]["dir"]),  # trimmed files and reports will be saved in this directory
        user_args=config["prep_args"]["trim_galore_se"],  # user args can be adjusted in the config file
        fastqc_rep_dir = expand("{rep_dir}/02_post_trim_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),  # fastqc reports will be moved to this directory
        trim_rep_dir = expand("{rep_dir}/02_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"])  # trim reports will be moved to this directory

    wildcard_constraints:
        layout="pe"

    shell:
        """
        mkdir -p {params.output_dir}
        trim_galore {params.user_args} --paired --fastqc --output_dir {params.output_dir} {input.r1} {input.r2} > {log} 2>&1
        mkdir -p {params.fastqc_rep_dir}
        mv -t {params.fastqc_rep_dir} {params.output_dir}/{wildcards.accession}_1_trimmed_fastqc.html {params.output_dir}/{wildcards.accession}_1_trimmed_fastqc.zip {params.output_dir}/{wildcards.accession}_2_trimmed_fastqc.html {params.output_dir}/{wildcards.accession}_2_trimmed_fastqc.zip
        mkdir -p {params.trim_rep_dir}
        mv -t {params.trim_rep_dir} {params.output_dir}/{wildcards.accession}.fastq_trimming_report.txt
        """