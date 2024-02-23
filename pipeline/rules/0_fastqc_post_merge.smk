### fastqc rules post merge ###
# added sleep 5 due to latency issues


rule fastqc_post_merge:
	input: 
		expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"])

	output:
		expand("{rep_dir}/fastqc/post-merge/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../pre-processing/logs/rule-logs/fastqc_post_merge/{ref}/fastqc_post_merge-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/fastqc/post-merge", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""



