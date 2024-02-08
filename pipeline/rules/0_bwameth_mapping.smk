#outputs not added to snakefile yet
#add inport of rule and expected outputs to snakefile for piped version to test
# outputs need to be created for se and pe versions of the piped rule
# STANDARD OUT WILL CAPTURE ALL OUTPUTS BECAUSE BWAMETH.PY SENDS OUTPUT FILE TO STDOUT

# add following outputs to snakefile before testing.
# expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{sample.ref}/{sample.patient_id}/{sample.group}-{sample.srx_id}/{sample.accession}{suf}", data_dir=config["data"]["dir"],sample=sample_info_se.itertuples(), suf=[".bam"])
# expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{sample.ref}/{sample.patient_id}/{sample.group}-{sample.srx_id}/{sample.accession}{suf}", data_dir=config["data"]["dir"],sample=sample_info_pe.itertuples(), suf=[".bam"])


rule bwameth_mapping_se:
    input:
        read_se = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa.gz", bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"])
            
    output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"]),
        bwa_report = expand("{rep_dir}/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"])
         
    log: 
        "../pre-processing/logs/rule-logs/bwameth_mapping_se/{ref}/bwameth_mapping_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
        
    conda: 
        "../env/bwameth.yaml"

    threads: 3

    params: 
        accession = "{accession}",
        sam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.sam", data_dir=config["data"]["dir"])

    shell:
        """
        echo "aligning {input.read_se} to {input.index}" > {log}
        bwameth.py --threads {threads} --reference {input.index} {input.read_se} > {params.sam} 2> {output.bwa_report}
        echo "converting {params.sam} to {output.bam}" >> {log}
        samtools view -S -b {params.sam} > {output.bam}
        echo "alignment and conversion complete for {params.accession}" >> {log} 
        echo "removing {params.sam} file to save space" >> {log}
        rm {params.sam} >> {log} 2>> {log}
        """

# # test logging
# rule bwameth_mapping_se_pipe:
#     input:
#         read_se = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", data_dir=config["data"]["dir"]),
#         index = expand("{bwa_idx_dir}/{fasta}.fa.gz", bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"])

#     output:
#         bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"]) 

#     log:
#         "../pre-processing/logs/rule-logs/bwameth_mapping_se_pipe/{ref}/bwameth_mapping_se_pipe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

#     conda:
#         "../env/bwameth.yaml"

#     threads: 6

#     params: 
#         out_dir = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data"]["dir"])

#     wildcard_constraints:
#         layout = "se"

#     shell:
#         """
#         mkdir -p {params.out_dir}
#         bwameth.py --threads {threads} --reference {input.index} {input.read_se} 2> {log} | samtools view -S -b > {output.bam}
#         """

rule bwameth_mapping_pe_pipe:
    input:
        read_1 = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_trimmed.fq", data_dir=config["data"]["dir"]),
        read_2 = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_trimmed.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa.gz", bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"])

    output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"]) 

    log:
         "../pre-processing/logs/rule-logs/bwameth_mapping_pe_pipe/{ref}/bwameth_mapping_se_pipe-{ref}-{patient_id}-{group}-{layout}-{srx_id}-{accession}.log"

    conda:
        "../env/bwameth.yaml"

    threads: 6

    params: 
        out_dir = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", data_dir=config["data"]["dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        mkdir -r {params.out_dir}
        bwameth.py --threads {threads} --reference {input.index} {input.read_1} {input.read_2} | samtools view -S -b > {output.bam}
        """
