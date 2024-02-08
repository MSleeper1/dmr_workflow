# trim_galore outputs fastqc reports into the trimming directory. This rule will simply move the files to the reports directory

rule move_reports_se:
    input:
        trim_report = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", data_dir=config["data"]["dir"]),
        fastqc = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", data_dir=config["data"]["dir"], suf = ["html", "zip"])

    output:
        trim_report = expand("{rep_dir}/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", rep_dir=config["reports_dir"]),
        fastqc = expand("{rep_dir}/fastqc/post-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], suf = ["html", "zip"])

    params:
        qc_rep_dir = expand("{rep_dir}/fastqc/post-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),
        trim_rep_dir = expand("{rep_dir}/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"])

    wildcard_constraints:
        layout = "se"

    shell:
        """
        mkdir -p {params.qc_rep_dir}
        mv -t {params.qc_rep_dir} {input.fastqc}
        mkdir -p {params.trim_rep_dir}
        mv {input.trim_report} {output.trim_report}
        """

rule move_reports_pe:
    input:
        trim_report = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", data_dir=config["data"]["dir"]),
        fastqc = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}{read}_trimmed_fastqc.{suf}", data_dir=config["data"]["dir"], read=["_1", "_2" ], suf = ["html", "zip"])

    output:
        trim_report = expand("{rep_dir}/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq_trimming_report.txt", rep_dir=config["reports_dir"]),
        fastqc = expand("{rep_dir}/fastqc/post-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}{read}_trimmed_fastqc.{suf}", rep_dir=config["data"]["dir"], read=["_1", "_2" ], suf = ["html", "zip"])

    params:
        qc_rep_dir = expand("{rep_dir}/fastqc/post-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),
        trim_rep_dir = expand("{rep_dir}/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        mkdir -p {params.qc_rep_dir}
        mv -t {params.qc_rep_dir} {input.fastqc}
        mkdir -p {params.trim_rep_dir}
        mv {input.trim_report} {output.trim_report}
        """