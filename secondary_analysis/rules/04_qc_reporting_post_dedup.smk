
### Snakemake rules for quality control reporting of trimmed, aligned, and deduplicated sequence files (03) ###
# fastqc and samtools stats are run on the deduplicated sequence files

### FASTQC RULES ###
# input: trimmed, alignned, and deduplicated sequence files (bam)
# output: fastqc reports for deduplicated sequence files (html and zip)

# Rule to run fastqc on bam files after deduplication
rule fastqc_post_dedup_bwa:
    input: 
        bwa_bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir = config["data_dir"]),
    
    output:
        expand("{root}/{rep_dir}/04_fastqc_post_dedup_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    log:
        "logs/secondary_rules/04_fastqc_post_dedup_bwa/04_fastqc_post_dedup-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    params:
        output_dir = expand("{root}/{rep_dir}/04_fastqc_post_dedup_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir = config["reports_dir"])

    shell: 
        """
		mkdir -p {params.output_dir}
		echo "Running fastqc on {input.bwa_bam}" > {log}
		fastqc -o {params.output_dir} {input.bwa_bam} >> {log} 2>&1
		echo "Done" >> {log}
		"""

rule fastqc_post_dedup_bis:
    input: 
        bis_bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"])

    output:
        expand("{root}/{rep_dir}/04_fastqc_post_dedup_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

    log:
        "logs/secondary_rules/04_fastqc_post_dedup_bis/04_fastqc_post_dedup-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/fastqc.yaml"

    params:
        output_dir = expand("{root}/{rep_dir}/04_fastqc_post_dedup_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir = config["reports_dir"])

    shell: 
        """
		mkdir -p {params.output_dir}
		echo "Running fastqc on {input.bis_bam}" > {log}
		fastqc -o {params.output_dir} {input.bis_bam} >> {log} 2>&1
		echo "Done" >> {log}
		"""

### SAMTOOLS STATS RULE ###
# samtools stats is a program that generates general statistics for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: samtools stats report for deduplicated sequence files (.stats text file)
rule samtools_stats_bwa:
    input:
        bam = expand("{root}/{data_dir}/04_deduped_sambamba/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", root = config["root"], data_dir=config["data_dir"]) # sambamba output

    output:
        report = expand("{root}/{rep_dir}/04_samtools_post_dedup_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.stats", root = config["root"], rep_dir=config["reports_dir"])

    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools stats on {input.bam}" > {output.report}
        samtools stats -p -d {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """

rule samtools_stats_bis:
    input:
        bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"])

    output:
        report = expand("{root}/{rep_dir}/04_samtools_post_dedup_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam.stats", root = config["root"], rep_dir=config["reports_dir"])

    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools stats on {input.bam}" > {output.report}
        samtools stats -p -d {input.bam} >> {output.report}
        echo "Done" >> {output.report}
        """


