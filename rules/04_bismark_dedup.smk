rule deduplicate_bismark:
    input: 
        expand("{data_dir}/03_aligned_bismark_bwt2/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2.bam", data_dir=config["data_dir"])
    
    output:
        bam = temporary(expand("{data_dir}/04_bismark_deduped/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_deduped_bismark_bt2.bam", data_dir=config["data_dir"])),
        report = expand("{rep_dir}/04_bismark_deduplication/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_bismark_deduplication_report.txt", rep_dir=config["reports_dir"])
    
    log:
        "../pre-processing/logs/rule-logs/04_deduplicate_bismark/{ref}/04_deduplicate_bismark-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
    
    conda:
        "../environment_files/bismark.yaml"

    shadow:
        "shallow"
    
    params:
        extra=""  # optional params string
    
    wrapper:
        "v3.4.1/bio/bismark/deduplicate_bismark"