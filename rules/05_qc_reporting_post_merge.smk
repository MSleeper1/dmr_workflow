### Snakemake rules for quality control of merged bam files (05) ###
# qualimap, samtools stats, fastqc, and featureCounts are used to generate quality control reports for the merged bam files.
# reports are saved to the reports directory. 
# directories will have 05 prefix to indicate that they are part of the 05_quality_control section of the pipeline reporting on merged bam files.

### QUALIMAP ###
# Qualimap is used to generate quality control reports for the merged bam files.
# Qualimaps flag --outdir is a relative path. To avoid issues with this, the ouput will be moved to the reports directory after qualimap is finished.

# qualimap_post_merge: qualimap is run on the merged bam files. The output is moved to the reports directory after qualimap is finished.
# input: merged bam files and gtf file (must be unzipped)
# output: qualimap reports

rule qualimap_post_merge_bwa:
    input:
        expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]),
        gtf = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])

    output:
        directory(expand("{rep_dir}/05_qualimap_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged", rep_dir=config["reports_dir"]))
        
    log:
        "../pre-processing/logs/rule-logs/05_qualimap_post_merge_bwa/{ref}/05_qualimap_post_merge_bwa-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../environment_files/qualimap.yaml"

    params:
        temp_out = "{ref}-{patient_id}-{group}-{srx_id}-{layout}_merged"

    shell:
        """
        qualimap bamqc -bam {input} -c -sd -os -gd hg38 -gff {input.gtf} --outdir {params.temp_out} > {log} 2>&1
        mv -f -v {params.temp_out} {output} >> {log} 2>&1
        rm -rf -v {params.temp_out} >> {log} 2>&1
        """

rule qualimap_post_merge_bis:
    input:
        expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]),
        gtf = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])

    output:
        directory(expand("{rep_dir}/05_qualimap_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged", rep_dir=config["reports_dir"]))
        
    log:
        "../pre-processing/logs/rule-logs/05_qualimap_post_merge_bis/{ref}/05_qualimap_post_merge_bis-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../environment_files/qualimap.yaml"

    params:
        temp_out = "bis_{ref}-{patient_id}-{group}-{srx_id}-{layout}_merged"

    shell:
        """
        qualimap bamqc -bam {input} -c -sd -os -gd hg38 -gff {input.gtf} --outdir {params.temp_out} > {log} 2>&1
        mv -f -v {params.temp_out} {output} >> {log} 2>&1
        rm -rf -v {params.temp_out} >> {log} 2>&1
        """


### FEATURE COUNTS ###
# input.sam can be a bam file but must be called input.sam to function with wrapper

rule feature_counts_bwa:
    input:
        sam = expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]),
        annotation = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    
    output:
        expand("{rep_dir}/05_feature_counts_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "../pre-processing/logs/rule-logs/05_feature_counts_bwa/{ref}/05_feature_counts_bwa-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    threads:
        3

    conda:
        "../environment_files/feature_counts.yaml"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule feature_counts_bis:
    input:
        sam = expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]),
        annotation = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    
    output:
        expand("{rep_dir}/05_feature_counts_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "../pre-processing/logs/rule-logs/05_feature_counts_bis/{ref}/05_feature_counts_bis-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    threads:
        3

    conda:
        "../environment_files/feature_counts.yaml"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"


### FASTQC ###
# added sleep 5 due to latency issues


rule fastqc_post_merge_bwa:
	input: 
		expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"])

	output:
		expand("{rep_dir}/05_fastqc_post_merge_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../pre-processing/logs/rule-logs/05_fastqc_post_merge_bwa/{ref}/05_fastqc_post_merge-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

	conda:
		"../environment_files/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/05_fastqc_post_merge", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""

rule fastqc_post_merge_bis:
	input: 
		expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"])

	output:
		expand("{rep_dir}/05_fastqc_post_merge_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../pre-processing/logs/rule-logs/05_fastqc_post_merge_bis/{ref}/05_fastqc_post_merge-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

	conda:
		"../environment_files/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/05_fastqc_post_merge", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""




### SAMTOOLS STATS ###

rule samtools_stats_post_merge_bwa:
    input:
        bam = expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"])
   
    output:
        report = expand("{rep_dir}/05_samtools_post_merge_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.stats", rep_dir=config["reports_dir"])

    conda:
        "../environment_files/samtools.yaml"

    shell:
        """
        samtools stats -p -d {input.bam} > {output.report}
        """

rule samtools_stats_post_merge_bis:
    input:
        bam = expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"])

    output:
        report = expand("{rep_dir}/05_samtools_post_merge_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.stats", rep_dir=config["reports_dir"])

    conda:
        "../environment_files/samtools.yaml"

    shell:
        """
        samtools stats -p -d {input.bam} > {output.report}
        """

### MOSDEPTH RULE ###
# mosdepth is a program that calculates the depth of coverage for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: mosdepth report for deduplicated sequence files (global distribution text file, per-base bed file, and summary text file)

rule mosdepth_bwa:
    input:
        bam = expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]),
        bai = expand("{data_dir}/05_merged_sambamba_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.bai", data_dir=config["data_dir"])
    output:
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.global.dist.txt", rep_dir=config["reports_dir"]),
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.per-base.bed.gz", rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        summary = expand("{rep_dir}/05_mosdepth_post_merge_bwa/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.summary.txt", rep_dir=config["reports_dir"]) # this named output is required for prefix parsing
    log:
        "../pre-processing/logs/rule-logs/05_mosdepth_bwa/{ref}/05_mosdepth_bis-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"
    params:
        extra="--fast-mode -Q 10",  # optional
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v3.4.1/bio/mosdepth"

rule mosdepth_bis:
    input:
        bam = expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data_dir"]),
        bai = expand("{data_dir}/05_merged_sambamba_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.bai", data_dir=config["data_dir"])
    output:
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.global.dist.txt", rep_dir=config["reports_dir"]),
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.per-base.bed.gz", rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        summary = expand("{rep_dir}/05_mosdepth_post_merge_bis/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.summary.txt", rep_dir=config["reports_dir"]) # this named output is required for prefix parsing
    log:
        "../pre-processing/logs/rule-logs/05_mosdepth_bis/{ref}/05_mosdepth_bis-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"
    params:
        extra="--fast-mode -Q 10",  # optional
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v3.4.1/bio/mosdepth"
