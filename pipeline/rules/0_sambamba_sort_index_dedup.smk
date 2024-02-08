

rule sambamba_sort_index_markdups:
    input: 
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"])

    output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam", data_dir=config["data"]["dir"]),
        bai = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.bai", data_dir=config["data"]["dir"]),
        in_bam_report = expand("{rep_dir}/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam.flagstat", rep_dir=config["reports_dir"]),
        out_bam_report = expand("{rep_dir}/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_sorted_dedup.bam.flagstat", rep_dir=config["reports_dir"])

    log:
        stdout = "../pre-processing/logs/rule-logs/sambamba_sort_index_markdups/{ref}/sambamba_sort_index_markdups-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.out",
        stderr = "../pre-processing/logs/rule-logs/sambamba_sort_index_markdups/{ref}/sambamba_sort_index_markdups-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"

    conda:
        "../env/sambamba.yaml"

    params:
        temp_dir = expand("{data_dir}/tmp/{{ref}}/{{accession}}", data_dir=config["data"]["dir"]),
        rep_dir = expand("{rep_dir}/sambamba/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["data"]["dir"])

    threads: 6

    shell:
        """
        sambamba flagstat -t {threads} {input.bam} > {output.in_bam_report}
        sambamba sort -p -t {threads} -o {output.bam} --tmpdir {params.temp_dir} {input.bam} >2 {log.stderr} | \
        sambamba markdup -p -t {threads} --overflow-list-size 610000000 --tmpdir {params.temp_dir} {output.bam} >> {log.stdout} 2>> {log.stderr}
        sambamba index -p -t {threads} {output.bam} {output.bai} >> {log.stdout} 2>> {log.stderr}
        sambamba flagstat -t {threads} {output.bam} > {output.out_bam_report}
        """
