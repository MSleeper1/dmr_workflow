rule bismark_methylation_extractor_se:
    input: 
        expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])
    
    output:
        mbias_r1 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], rep_dir=config["reports_dir"]),
        # Only for PE BAMS:
        # mbias_r2="qc/meth/{sample}.M-bias_R2.png",

        mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
        splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"]),

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylone_CpG_cov = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], rep_dir=config["reports_dir"]),
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], rep_dir=config["reports_dir"]),

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CpG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),

        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        read_base_meth_state_chh = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHH_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"])
    
    log:
        "logs/secondary_rules/06_bismark_methylation_extractor_se/06_bismark_methylation_extractor_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/bismark.yaml"
    # shadow:
    #     "shallow"
    
    wildcard_constraints:
        layout="se"
    
    params:
        output_dir="meth",  # optional output dir
        extra="--gzip --comprehensive --bedGraph"  # optional params string
    
    wrapper:
        "v3.4.1/bio/bismark/bismark_methylation_extractor"


rule bismark_methylation_extractor_pe:
    input: 
        expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])
    
    output:
        mbias_r1 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], rep_dir=config["reports_dir"]),
        # Only for PE BAMS:
        mbias_r2 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R2.png", root = config["root"], rep_dir=config["reports_dir"]),

        mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
        splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"]),

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylone_CpG_cov = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], rep_dir=config["reports_dir"]),
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], rep_dir=config["reports_dir"]),

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CpG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),

        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHG_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"]),
        read_base_meth_state_chh = expand("{root}/{data_dir}/06_bismark_methyl_extractor/CHH_context_{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.txt.gz", root = config["root"], data_dir=config["data_dir"])
    
    log:
        "logs/secondary_rules/06_bismark_methylation_extractor_pe/06_bismark_methylation_extractor_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/bismark.yaml"
    # shadow:
    #     "shallow"    
    
    wildcard_constraints:
        layout="pe"
    
    params:
        output_dir="meth",  # optional output dir
        extra="--gzip --comprehensive --bedGraph"  # optional params string
    
    wrapper:
        "v3.4.1/bio/bismark/bismark_methylation_extractor"