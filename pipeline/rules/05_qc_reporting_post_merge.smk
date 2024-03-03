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

rule qualimap_post_merge:
    input:
        expand("{data_dir}/05_merged_sambamba/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"]),
        gtf = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])

    output:
        directory(expand("{rep_dir}/05_qualimap/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged", rep_dir=config["reports_dir"]))
        
    log:
        "../pre-processing/logs/rule-logs/05_qualimap_post_merge/{ref}/05_qualimap_post_merge-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../env/qualimap.yaml"

    params:
        temp_out = "{ref}-{patient_id}-{group}-{srx_id}-{layout}_merged"

    shell:
        """
        qualimap bamqc -bam {input} -c -sd -os -gd hg38 -gff {input.gtf} --outdir {params.temp_out} > {log} 2>&1
        mv -f -v {params.temp_out} {output} >> {log} 2>&1
        rm -rf -v {params.temp_out} >> {log} 2>&1
        """


### FEATURE COUNTS ###
# input.sam can be a bam file but must be called input.sam to function with wrapper

rule feature_counts:
    input:
        sam = expand("{data_dir}/05_merged_sambamba/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"]),
        annotation = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    
    output:
        expand("{rep_dir}/05_feature_counts/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "../pre-processing/logs/rule-logs/05_feature_counts/{ref}/05_feature_counts-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    threads:
        3

    conda:
        "../env/feature_counts.yaml"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"


### FASTQC ###


# added sleep 5 due to latency issues


rule fastqc_post_merge:
	input: 
		expand("{data_dir}/05_merged_sambamba/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"])

	output:
		expand("{rep_dir}/05_fastqc_post_merge/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"../pre-processing/logs/rule-logs/05_fastqc_post_merge/{ref}/05_fastqc_post_merge-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

	conda:
		"../env/fastqc.yaml"

	params:
		output_dir = expand("{rep_dir}/05_fastqc_post_merge", rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input} > {log} 2>&1
		"""




### SAMTOOLS STATS ###

rule samtools_stats_post_merge:
    input:
        bam = expand("{data_dir}/05_merged_sambamba/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"])

    output:
        report = expand("{rep_dir}/05_samtools_post_merge/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.stats", rep_dir=config["reports_dir"])

    conda:
        "../env/samtools.yaml"

    shell:
        """
        samtools stats -p -d {input.bam} > {output.report}
        """