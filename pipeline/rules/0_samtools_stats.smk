
rule samtools_stats:
    input:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", data_dir=config["data"]["dir"]) # sambamba output

    output:
        report = expand("{rep_dir}/samtools/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_dedup.stats", rep_dir=config["reports_dir"])

    conda:
        "../env/samtools.yaml"

    shell:
        """
        samtools stats -p -d {input.bam} > {output.report}
        """
