### BWA-METH alignment rules ###
# bwa-meth is a tool for aligning bisulfite-treated reads to a reference genome. It is a wrapper around BWA and is designed to work with WGBS data.
# samtools is used to convert the sam file to a bam file and the sam file is removed to save space
# input: trimmed fastq files and reference genome
# output: aligned bam files and alignment report (captured standard error report of bwameth)
# LOGGING STANDARD OUT OF BWA-METH WILL CAPTURE ALL OUTPUTS BECAUSE BWAMETH.PY SENDS OUTPUT FILE TO STDOUT

# bwameth_mapping_se_pipe: align single-end reads to reference genome using bwameth
rule bwameth_mapping_se:
    input:
        read_se = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa.gz", bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"])
            
    output:
        bam = temporary(expand("{data_dir}/03_aligned_bwameth/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"])),
        bwa_report = expand("{rep_dir}/03_bwameth/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"])
         
    log: 
        "../pre-processing/logs/rule-logs/03_bwameth_mapping_se/{ref}/03_bwameth_mapping_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
        
    conda: 
        "../env/bwameth.yaml"

    threads: 3

    params: 
        accession = "{accession}",
        sam = expand("{data_dir}/03_aligned_bwameth/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.sam", data_dir=config["data"]["dir"])

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

# bwameth_mapping_pe_pipe: align paired-end reads to reference genome using bwameth
rule bwameth_mapping_pe_pipe:
    input:
        read_1 = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_trimmed.fq", data_dir=config["data"]["dir"]),
        read_2 = expand("{data_dir}/02_trimmed_trim_galore/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_trimmed.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa.gz", bwa_idx_dir=config["ref"]["bwa_idx_dir"], fasta=config["ref"]["fasta"])

    output:
        bam = temporary(expand("{data_dir}/03_aligned_bwameth/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.bam", data_dir=config["data"]["dir"])),
        bwa_report = expand("{rep_dir}/03_bwameth/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"])
         
    log:
        "../pre-processing/logs/rule-logs/03_bwameth_mapping_pe_pipe/{ref}/03_bwameth_mapping_se_pipe-{ref}-{patient_id}-{group}-{layout}-{srx_id}-{accession}.log"

    conda:
        "../env/bwameth.yaml"

    threads: 6

    params: 
        accession = "{accession}",
        sam = expand("{data_dir}/03_aligned_bwameth/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.sam", data_dir=config["data"]["dir"])

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        echo "aligning {input.read_1} and {input.read_2} to {input.index}" > {log}
        bwameth.py --threads {threads} --reference {input.index} {input.read_1} {input.read_2} > {params.sam} 2> {output.bwa_report}
        echo "converting {params.sam} to {output.bam}" >> {log}
        samtools view -S -b {params.sam} > {output.bam}
        echo "alignment and conversion complete for {params.accession}" >> {log} 
        echo "removing {params.sam} file to save space" >> {log}
        rm {params.sam} >> {log} 2>> {log}
        """
