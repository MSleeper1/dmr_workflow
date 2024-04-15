rule bismark_deduplicate:
    input: 
        expand("{data_dir}/03_aligned_bismark_bwt2/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2.bam", data_dir=config["data_dir"])
    
    output:
        bam = temporary(expand("{data_dir}/04_bismark_deduped/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", data_dir=config["data_dir"])),
        report = expand("{rep_dir}/04_bismark_deduplication/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplication_report.txt", rep_dir=config["reports_dir"])
    
    log:
        "../pre-processing/logs/rule-logs/04_bismark_deduplicate/{ref}/04_bismark_deduplicate-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
    
    conda:
        "../environment_files/bismark.yaml"

    shadow:
        "shallow"
    
    params:
        output_dir = expand("{data_dir}/04_bismark_deduped/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data_dir"]),
        base_name = "{accession}_bismark",
        report_dir = expand("{rep_dir}/04_bismark_deduplication/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"]),
        report = "{accession}_bismark.deduplication_report.txt"

    shell:
        """
        deduplicate_bismark --bam {input} --output_dir {params.output_dir} --outfile {params.base_name} > {log} 2> {log}
        mv -f -v --target-directory={params.report_dir} {params.output_dir}/{params.report} >> {log} 2>> {log}
        """