rule move_bis_reports_se:
    input:
        mbias_r1 = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], data_dir=config["data_dir"]),
        mbias_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
        splitting_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),
        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylone_CpG_cov = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], data_dir=config["data_dir"])
    
    output:
        mbias_r1 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], rep_dir=config["reports_dir"]),
        mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
        splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        methylone_CpG_cov = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], rep_dir=config["reports_dir"]),
        methylome_CpG_mlevel_bedGraph = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], rep_dir=config["reports_dir"])
    
    log:
        "logs/secondary_rules/06_move_bis_reports/06_move_bis_reports_se-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    # shadow:
    #     "shallow"    
    
    wildcard_constraints:
        layout="se"
    
    params:
        output_dir = expand("{root}/{rep_dir}/06_bismark_methyl_extractor", root = config["root"], rep_dir=config["reports_dir"])

    shell:
        """
        mkdir -p {params.output_dir}
        mv {input.mbias_r1} {output.mbias_r1}
        mv {input.mbias_report} {output.mbias_report}
        mv {input.splitting_report} {output.splitting_report}
        mv {input.methylone_CpG_cov} {output.methylone_CpG_cov}
        mv {input.methylome_CpG_mlevel_bedGraph} {output.methylome_CpG_mlevel_bedGraph}
        """

    
rule move_bis_reports_pe:
    input:
        mbias_r1 = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], data_dir=config["data_dir"]),
        # Only for PE BAMS:
        mbias_r2 = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R2.png", root = config["root"], data_dir=config["data_dir"]),
        mbias_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], data_dir=config["data_dir"]),
        splitting_report = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], data_dir=config["data_dir"]),
        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylone_CpG_cov = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], data_dir=config["data_dir"]),
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph = expand("{root}/{data_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], data_dir=config["data_dir"])
    
    output:
        mbias_r1 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R1.png", root = config["root"], rep_dir=config["reports_dir"]),
        mbias_r2 = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias_R2.png", root = config["root"], rep_dir=config["reports_dir"]),
        mbias_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.M-bias.txt", root = config["root"], rep_dir=config["reports_dir"]),
        splitting_report = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_splitting_report.txt", root = config["root"], rep_dir=config["reports_dir"]),
        methylone_CpG_cov = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bismark.cov.gz", root = config["root"], rep_dir=config["reports_dir"]),
        methylome_CpG_mlevel_bedGraph = expand("{root}/{rep_dir}/06_bismark_methyl_extractor/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}.bedGraph.gz", root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/06_move_bis_reports/06_move_bis_reports_pe-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    wildcard_constraints:
        layout="pe"

    params:
        output_dir = expand("{root}/{rep_dir}/06_bismark_methyl_extractor", root = config["root"], rep_dir=config["reports_dir"])
    
    shell:
        """
        mkdir -p {params.output_dir}
        mv {input.mbias_r1} {output.mbias_r1}
        mv {input.mbias_r2} {output.mbias_r2}
        mv {input.mbias_report} {output.mbias_report}
        mv {input.splitting_report} {output.splitting_report}
        mv {input.methylone_CpG_cov} {output.methylone_CpG_cov}
        mv {input.methylome_CpG_mlevel_bedGraph} {output.methylome_CpG_mlevel_bedGraph}
        """

