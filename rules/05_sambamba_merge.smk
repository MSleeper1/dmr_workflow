
### sambamba merge rule for deduped bam files ###
# sambamba is a tool for working with SAM/BAM files. It is faster than samtools and can be used to merge bam files.

# sambamba_merge: determines if merging is needed based on SRX ID. 
# If only one file is associated with the SRX ID, a symlink is created. 
# If more than one file is associated with the SRX ID, the files are merged using sambamba merge. 
# The merged file is then used for the next step in the pipeline.
# input: deduped bam files (04)
# output: merged bam file (05)

# # has been tested with files that do not need to be merged
# # still need to test rule for sambamba merge with files that need to be merged
# symlinks were experiencing latency wait time errors. To avoid specifying latency wait time for all jobs, I added a sleep 10 seconds command to the shell script. This will allow the symlink to be created before the next rule_all is executed.

rule sambamba_merge_bwameth:
    input:
       bams = lambda wildcards: expand("{data_dir}/04_deduped_sambamba/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{accession}_trimmed_sorted_dedup.bam", data_dir=config["data_dir"], accession = sample_info[sample_info["srx_id"] == wildcards.srx_id]["accession"].tolist())
    
    output:
        merged_bam = temporary(expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"])),
        bai = expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.bai", data_dir=config["data_dir"])

    log:
        "../pre-processing/logs/rule-logs/05_sambamba_merge_bwameth/{ref}/05_sambamba_merge_bwameth-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../environment_files/sambamba.yaml"

    shadow:
        "shallow"

    threads: 3
    
    wildcard_constraints:
        srx_id = "|".join(sample_info["srx_id"].tolist()),
        accession = "|".join(sample_info["accession"].tolist())

    shell:
        """
        FILES=()
        for i in {input.bams}; do
            FILES+=($i)
        done
        input_len=${{#FILES[@]}}
        if [ $input_len -gt 1 ]; then
            echo "Merging files: {input.bams} into {output.merged_bam} with sambamba merge..." > {log}
            sambamba merge -t {threads} {output.merged_bam} {input.bams} >> {log} 2>> {log}
            sambamba index -t {threads} {output.merged_bam} {output.bai} >> {log} 2>> {log}
        else
            echo "Only one file, no need to merge. Creating symlink for file at expected output location: {output.merged_bam}" > {log}
            cp -f {input.bams} {output.merged_bam}
            cp -f {input.bams}.bai {output.bai}
        fi
        sleep 10
        """

rule sambamba_merge_bismark:
    input:
        bams = lambda wildcards: expand("{data_dir}/04_bismark_deduped/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{accession}_bismark.deduplicated.bam", data_dir=config["data_dir"], accession = sample_info[sample_info["srx_id"] == wildcards.srx_id]["accession"].tolist())
    
    output:
        merged_bam = temporary(expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"])),
        bai = expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.bai", data_dir=config["data_dir"])

    log:
        "../pre-processing/logs/rule-logs/05_sambamba_merge_bismark/{ref}/05_sambamba_merge_bismark-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../environment_files/sambamba.yaml"

    shadow:
        "shallow"

    threads: 3
    
    wildcard_constraints:
        srx_id = "|".join(sample_info["srx_id"].tolist()),
        accession = "|".join(sample_info["accession"].tolist())

    shell:
        """
        FILES=()
        for i in {input.bams}; do
            FILES+=($i)
        done
        input_len=${{#FILES[@]}}
        if [ $input_len -gt 1 ]; then
            echo "Merging files: {input.bams} into {output.merged_bam} with sambamba merge..." > {log}
            sambamba merge -t {threads} {output.merged_bam} {input.bams} >> {log} 2>> {log}
            sambamba index -t {threads} {output.merged_bam} {output.bai} >> {log} 2>> {log}
        else
            echo "Only one file, no need to merge. Creating symlink for file at expected output location: {output.merged_bam}" > {log}
            cp -f {input.bams} {output.merged_bam}
            cp -f {input.bams}.bai {output.bai}
        fi
        sleep 10
        """