
### Snakemake rules for fastqc quality control of raw sequence files ###

# Rule to run fastqc on single-end sequence files
rule fastqc_se:
	input: 
		expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.fastq", data_dir = config["data"]["dir"])

	output:
		expand("{rep_dir}/quality/fastqc/pretrim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../logs/fastqc_se/fastqc_se-{ref}-{patient_id}-{group}-{srx_id}-{accession}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/quality/fastqc/pretrim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""

# Rule to run fastqc on paired-end sequence files
rule fastqc_pe:
	input: 
		r1 = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_1.fastq", data_dir = config["data"]["dir"]),
		r2 = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_2.fastq", data_dir = config["data"]["dir"])

	output:
		r1 = expand("{rep_dir}/quality/fastqc/pretrim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_1_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"]),
		r2 = expand("{rep_dir}/quality/fastqc/pretrim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_2_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
        "../logs/fastqc_pe/fastqc_pe-{ref}-{patient_id}-{group}-{srx_id}-{accession}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/quality/fastqc/pretrim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input.r1} {input.r2} > {log} 2>&1
		"""



