
### Snakemake rules for fastqc quality control of raw sequence files ###
# added sleep 5 due to latency issues

# Rule to run fastqc on single-end sequence files
rule fastqc_post_dedup:
	input: 
		expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", data_dir = config["data"]["dir"])

	output:
		expand("{rep_dir}/fastqc/post-dedup/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../pre-processing/logs/rule-logs/fastqc_se_post_dedup/{ref}/fastqc_se_post_dedup-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/fastqc/post-dedup/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""




