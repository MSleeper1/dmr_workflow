rule bismark_deduplicate:
    input: 
        expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2.bam", root = config["root"], data_dir=config["data_dir"]),
        
    output:
        bam = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplicated.bam", root = config["root"], data_dir=config["data_dir"]),
        report = expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark.deduplication_report.txt", root = config["root"], rep_dir=config["reports_dir"])
    
    log:
        "logs/secondary_rules/04_bismark_deduplicate/04_bismark_deduplicate-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
    
    conda:
        "../../environment_files/bismark.yaml"

    # shadow:
    #     "shallow"
    
    params:
        output_dir = expand("{root}/{data_dir}/04_bismark_deduped/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),
        base_name = "{accession}_bismark",
        report_dir = expand("{root}/{rep_dir}/04_bismark_deduplication/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        report = "{accession}_bismark.deduplication_report.txt"

    shell:
        """
        echo "making output directory {params.output_dir}" > {log}
        mkdir -p {params.output_dir} 2>>{log}
        echo "running bismark deduplication on {input}"
        deduplicate_bismark --bam {input} --output_dir {params.output_dir} --outfile {params.base_name} > {log} 2> {log}
        echo "done with deduplication"
        echo "making report directory: {params.report_dir}" >> {log}
        mkdir -p {params.report_dir} 2>>{log}
        echo "moving {params.output_dir}/{params.base_name}.deduplication_report.txt to {params.report_dir}" >> {log}
        mv -f -v --target-directory={params.report_dir} {params.output_dir}/{params.report} >> {log} 2>> {log}
        echo "done"
        """