
### Snakemake rules for fastqc quality control of raw sequence files ###

# Rule to run fastqc on single-end sequence files
rule fastqc_se:
	input: 
		expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir = config["data"]["dir"])

	output:
		expand("{rep_dir}/fastqc/pre-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../pre-processing/logs/rule-logs/fastqc_se/{ref}/fastqc_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/fastqc/pre-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir = config["reports_dir"])

	wildcard_constraints:
		layout = "se"

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""

# Rule to run fastqc on paired-end sequence files
rule fastqc_pe:
	input: 
		r1 = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", data_dir = config["data"]["dir"]),
		r2 = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", data_dir = config["data"]["dir"])

	output:
		r1 = expand("{rep_dir}/fastqc/pre-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"]),
		r2 = expand("{rep_dir}/fastqc/pre-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
        "../pre-processing/logs/rule-logs/fastqc_pe/{ref}/fastqc_pe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/fastqc/pre-trim/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir = config["reports_dir"])

	wildcard_constraints:
		layout = "pe"

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input.r1} {input.r2} > {log} 2>&1
		"""



