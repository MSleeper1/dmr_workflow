# flagstat outputs
        # in_bam_report = expand("{rep_dir}/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam.flagstat", rep_dir=config["reports_dir"]),
        # out_bam_report = expand("{rep_dir}/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.flagstat", rep_dir=config["reports_dir"])
        # rep_dir = expand("{rep_dir}/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["data"]["dir"])

rule sambamba_sort_index_markdups:
    input: 
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"])

    output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", data_dir=config["data"]["dir"]),
        bai = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.bai", data_dir=config["data"]["dir"]),
        log = expand("{rep_dir}/sambamba_sort_index_markdups/{{ref}}/sambamba_sort_index_markdups-{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}-{{accession}}.log", rep_dir=config["reports_dir"])

    log:
        expand("{rep_dir}/sambamba_sort_index_markdups/{{ref}}/sambamba_sort_index_markdups-{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}-{{accession}}.log", rep_dir=config["reports_dir"])

    conda:
        "../env/sambamba.yaml"

    params:
        temp_dir = expand("{data_dir}/temp/sambamba/{{ref}}/{{accession}}", data_dir=config["data"]["dir"]),
        sorted_bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted.bam", data_dir=config["data"]["dir"]),

    threads: 3

    shell:
        """
        mkdir -p {params.temp_dir}
        echo "Sorting bam file..." > {log}
        sambamba sort -t {threads} -o {params.sorted_bam} --tmpdir {params.temp_dir} {input.bam} >> {log} 2>> {log}
        echo "Sorting complete. Now marking duplicates..." >> {log} 
        sambamba markdup -t {threads} --remove-duplicates --tmpdir {params.temp_dir} {params.sorted_bam} {output.bam} >> {log} 2>> {log}
        echo "Marking duplicates complete. Now indexing..." >> {log}
        sambamba index -t {threads} {output.bam} {output.bai} >> {log} 2>> {log}
        echo "Indexing complete. Now removing intermediate bam file.." >> {log}
        rm -f {params.sorted_bam} >> {log} 2>> {log}
        echo "Intermediate bam file removed. Done." >> {log}
        """
