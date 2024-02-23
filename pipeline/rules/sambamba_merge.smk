rule sambamba_merge:
    input:
        expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{accession}_trimmed_sorted_dedup.bam", data_dir=config["data"]["dir"], accession = sample_info[sample_info['srx_id']=={srx_id}]) # sambamba_dedup output
    
    output:
        expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}_trimmed_bwameth_deduped_merged.bam", data_dir=config["data"]["dir"]) # sambamba_merge output

    log:
        "logs/sambamba-merge/{sample}.log"

    threads: 1

    wrapper:
        sambamba merge {output} {input}