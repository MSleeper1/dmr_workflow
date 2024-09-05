### Bismark mapping rules ###
# Bismark is a tool for aligning bisulfite treated sequencing reads to a reference genome.
# input: trimmed reads and reference genome
# output: aligned bam file and report

# bismark mapping for single end reads
rule bismark_mapping_se:
    input:
        index = expand("{root}/{genomes_dir}/{genome}/bismark/Bisulfite_Genome/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        read_se = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed.fq", root = config["root"], data_dir=config["data_dir"])
    
    output:
        bam = expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2.bam", root = config["root"], data_dir=config["data_dir"]),
        report = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_SE_report.txt", root = config["root"], rep_dir=config["reports_dir"])
    
    log:
        "logs/secondary_rules/03_bismark_mapping_se/03_bismark_mapping_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../../environment_files/bismark.yaml"

    # shadow: 
    #     "shallow"

    params:
        output_dir = expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),
        bismark_idx_dir = expand("{root}/{genomes_dir}/{genome}/bismark", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        temp_dir = expand("{root}/{data_dir}/temp/bismark/{{ref}}--{{accession}}", root = config["root"], data_dir=config["data_dir"]),
        report_dir = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        accession = "{accession}"

    wildcard_constraints:
        layout = "se"

    shell:
        """
        echo "making output directory {params.output_dir}" > {log}
        mkdir -p {params.output_dir} 2>>{log}
        echo "making temp directory {params.temp_dir}" >> {log}
        mkdir -p {params.temp_dir} 2>>{log}
        echo "running bismark on {input.read_se}"
        bismark --bowtie2 --temp_dir {params.temp_dir} --output_dir {params.output_dir} {params.bismark_idx_dir} {input.read_se} >> {log} 2>>{log}
        echo "done with mapping"
        echo "making report directory: {params.report_dir}" >> {log}
        mkdir -p {params.report_dir} 2>>{log}
        echo "moving {params.output_dir}/{params.accession}_trimmed_bismark_bt2_SE_report.txt to {params.report_dir}" >> {log}
        mv -f -v --target-directory={params.report_dir} {params.output_dir}/{params.accession}_trimmed_bismark_bt2_SE_report.txt 2>>{log}
        rm -rf {params.temp_dir}
        echo "done"
        """

# bismark mapping for paired end reads
rule bismark_mapping_pe:
    input:
        index = expand("{root}/{genomes_dir}/{genome}/bismark/Bisulfite_Genome/", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        r1 = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_trimmed.fq", root = config["root"], data_dir=config["data_dir"]),
        r2 = expand("{root}/{data_dir}/02_trimmed_trim_galore/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_trimmed.fq", root = config["root"], data_dir=config["data_dir"])
   
    output:
        bam = expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2.bam", root = config["root"], data_dir=config["data_dir"]),
        report = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_trimmed_bismark_bt2_PE_report.txt", root = config["root"], rep_dir=config["reports_dir"])
    
    log:
        "../pre-processing/logs/rule-logs/03_bismark_mapping_pe/{ref}/03_bismark_mapping_pe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"

    conda:
        "../environment_files/bismark.yaml"

    # shadow: 
    #     "shallow"

    params:
        output_dir = expand("{root}/{data_dir}/03_aligned_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], data_dir=config["data_dir"]),
        bismark_idx_dir = expand("{root}/{genomes_dir}/{genome}/bismark", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"]),
        temp_dir = expand("{root}/{data_dir}/temp/bismark/{{ref}}--{{accession}}", root = config["root"], data_dir=config["data_dir"]),
        report_dir = expand("{root}/{rep_dir}/03_bismark_bwt2/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", root = config["root"], rep_dir=config["reports_dir"]),
        accession = "{accession}"

    wildcard_constraints:
        layout = "pe"

    shell:
        """
        echo "making output directory {params.output_dir}" > {log}
        mkdir -p {params.output_dir} 2>>{log}
        echo "making temp directory {params.temp_dir}" >> {log}
        mkdir -p {params.temp_dir} 2>>{log}
        echo "running bismark on {input.r1} and {input.r2}"
        bismark --bowtie2 --temp_dir {params.temp_dir} --output_dir {params.output_dir} {params.bismark_idx_dir} -1 {input.r1} -2 {input.r2} >> {log} 2>>{log}
        echo "done with mapping"
        echo "making report directory: {params.report_dir}" >> {log}
        mkdir -p {params.report_dir} 2>>{log}
        echo "moving {params.output_dir}/{params.accession}_trimmed_bismark_bt2_SE_report.txt to {params.report_dir}" >> {log}
        mv -f -v --target-directory={params.report_dir} {params.output_dir}/{params.accession}_trimmed_bismark_bt2_PE_report.txt 2>>{log}
        rm -rf {params.temp_dir}
        echo "done"
        """