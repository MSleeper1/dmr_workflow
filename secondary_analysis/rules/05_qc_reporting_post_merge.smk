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
        expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
    
    output:
        directory(expand("{root}/{rep_dir}/05_qualimap_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged", root = config["root"], rep_dir=config["reports_dir"]))
        
    log:
        "logs/secondary_rules/05_qualimap_post_merge_bwa/05_qualimap_post_merge_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/qualimap.yaml"

    params:
        temp_out = "{ref}--{patient_id}-{group}-{srx_id}-{layout}_merged"

    shell:
        """
        echo "Running qualimap on {input}" > {log}
        qualimap bamqc -bam {input} -c -sd -os -gd hg38 -gff {input.gtf} --outdir {params.temp_out} >> {log} 2>&1
        echo "Moving qualimap output to reports directory" >> {log}
        mv -f -v {params.temp_out} {output} >> {log} 2>&1
        echo "removing temp output directory" >> {log}
        rm -rf -v {params.temp_out} >> {log} 2>&1
        echo "Done" >> {log}
        """

rule qualimap_post_merge_bis:
    input:
        expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        gtf = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
    
    output:
        directory(expand("{root}/{rep_dir}/05_qualimap_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged", root = config["root"], rep_dir=config["reports_dir"]))
        
    log:
        "logs/secondary_rules/05_qualimap_post_merge_bis/05_qualimap_post_merge_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../../environment_files/qualimap.yaml"

    params:
        temp_out = "bis_{ref}--{patient_id}-{group}-{srx_id}-{layout}_merged"

    shell:
        """
        echo "Running qualimap on {input}" > {log}
        qualimap bamqc -bam {input} -c -sd -os -gd hg38 -gff {input.gtf} --outdir {params.temp_out} >> {log} 2>&1
        echo "Moving qualimap output to reports directory" >> {log}
        mv -f -v {params.temp_out} {output} >> {log} 2>&1
        echo "removing temp output directory" >> {log}
        rm -rf -v {params.temp_out} >> {log} 2>&1
        echo "Done" >> {log}
        """


### FEATURE COUNTS ###
# input.sam can be a bam file but must be called input.sam to function with wrapper

rule feature_counts_bwa:
    input:
        sam = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/05_feature_counts_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/05_feature_counts_bwa/05_feature_counts_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads:
        3

    conda:
        "../../environment_files/feature_counts.yaml"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule feature_counts_bis:
    input:
        sam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        annotation = expand("{root}/{genomes_dir}/{genome}/{gtf}.gtf", root = config["root"], genomes_dir = config["genomes_dir"], genome = config["ref"]["genome"], gtf = config["ref"]["gtf"])
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    
    output:
        expand("{root}/{rep_dir}/05_feature_counts_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.featureCounts{suf}", root = config["root"], rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "logs/secondary_rules/05_feature_counts_bis/05_feature_counts_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

    threads:
        3

    conda:
        "../../environment_files/feature_counts.yaml"

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
		expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])

	output:
		expand("{root}/{rep_dir}/05_fastqc_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"logs/secondary_rules/05_fastqc_post_merge_bwa/05_fastqc_post_merge-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

	conda:
		"../../environment_files/fastqc.yaml"

	params:
		output_dir = expand("{root}/{rep_dir}/05_fastqc_post_merge", root = config["root"], rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
        echo "Running fastqc on {input}" > {log}
		fastqc -o {params.output_dir} {input} >> {log} 2>&1
        echo "Done" >> {log}
		"""

rule fastqc_post_merge_bis:
	input: 
		expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])

	output:
		expand("{root}/{rep_dir}/05_fastqc_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged_fastqc.{suf}", root = config["root"], rep_dir = config["reports_dir"], suf=["html","zip"])

	log:
		"logs/secondary_rules/05_fastqc_post_merge_bis/05_fastqc_post_merge-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"

	conda:
		"../../environment_files/fastqc.yaml"

	params:
		output_dir = expand("{root}/{rep_dir}/05_fastqc_post_merge", root = config["root"], rep_dir = config["reports_dir"])

	shell: 
		"""
		mkdir -p {params.output_dir}
        echo "Running fastqc on {input}" > {log}
		fastqc -o {params.output_dir} {input} >> {log} 2>&1
        echo "Done" >> {log}
		"""


### SAMTOOLS STATS ###

rule samtools_stats_post_merge_bwa:
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])
   
    output:
        report = expand("{root}/{rep_dir}/05_samtools_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.stats", root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/05_samtools_post_merge_bwa/05_samtools_post_merge_bwa-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools stats on {input}" > {log}
        samtools stats -p -d {input.bam} >> {output.report}
        echo "Done" >> {log}
        """

rule samtools_stats_post_merge_bis:
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"])

    output:
        report = expand("{root}/{rep_dir}/05_samtools_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.stats",  root = config["root"], rep_dir=config["reports_dir"])

    log:
        "logs/secondary_rules/05_samtools_post_merge_bis/05_samtools_post_merge_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    conda:
        "../../environment_files/samtools.yaml"

    shell:
        """
        echo "Running samtools stats on {input}" > {log}
        samtools stats -p -d {input.bam} >> {output.report}
        echo "Done" >> {log}
        """

### MOSDEPTH RULE ###
# mosdepth is a program that calculates the depth of coverage for sequence files
# input: trimmed, aligned, and deduplicated sequence files (bam)
# output: mosdepth report for deduplicated sequence files (global distribution text file, per-base bed file, and summary text file)

rule mosdepth_bwa:
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/05_merged_sambamba_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.bai", root = config["root"], data_dir=config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/05_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/05_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        summary = expand("{root}/{rep_dir}/05_mosdepth_post_merge_bwa/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing
    
    log:
        "logs/secondary_rules/05_mosdepth_bwa/05_mosdepth_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    params:
        extra="--fast-mode -Q 10",  # optional
    # additional decompression threads through `--threads`
    
    threads: 4  # This value - 1 will be sent to `--threads`
    
    wrapper:
        "v3.4.1/bio/mosdepth"

rule mosdepth_bis:
    input:
        bam = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", root = config["root"], data_dir=config["data_dir"]),
        bai = expand("{root}/{data_dir}/05_merged_sambamba_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam.bai", root = config["root"], data_dir=config["data_dir"])
    
    output:
        expand("{root}/{rep_dir}/05_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.global.dist.txt", root = config["root"], rep_dir=config["reports_dir"]),
        expand("{root}/{rep_dir}/05_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.per-base.bed.gz", root = config["root"], rep_dir=config["reports_dir"]), # produced unless --no-per-base specified
        summary = expand("{root}/{rep_dir}/05_mosdepth_post_merge_bis/{{ref}}--{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.mosdepth.summary.txt", root = config["root"], rep_dir=config["reports_dir"]) # this named output is required for prefix parsing
    
    log:
        "logs/secondary_rules/05_mosdepth_bis/05_mosdepth_bis-{ref}--{patient_id}-{group}-{srx_id}-{layout}.log"
    
    params:
        extra="--fast-mode -Q 10",  # optional
    # additional decompression threads through `--threads`
    
    threads: 4  # This value - 1 will be sent to `--threads`
    
    wrapper:
        "v3.4.1/bio/mosdepth"

# # Get the depth for each sample
# rule mosdepth:
#     input:
#         '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai',
#         bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
#     output:
#         '6_mosdepth/{sample}.sorted.markdupes.mosdepth.global.dist.txt',
#         '6_mosdepth/{sample}.sorted.markdupes.mosdepth.summary.txt',
#         '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz',
#         '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz.csi'
#     threads:
#         config['mosdepth']['threads']
#     params:
#         mapping_quality = config['mosdepth']['mapping_quality'],
#         mosdepth_path = config['paths']['mosdepth_path'],
#         out_prefix = '6_mosdepth/{sample}.sorted.markdupes'
#     shell:
#         '''
#         {params.mosdepth_path} \
#         -x \
#         -t {threads} \
#         -Q {params.mapping_quality} \
#         {params.out_prefix} \
#         {input.bam}
#         '''

# # Calculate the coverage from the mosdepth output
# rule calc_coverage:
#     input:
#         bed = '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz'
#     output:
#         '6_mosdepth/{sample}.sorted.markdupes.coverage.txt'
#     params:
#         genome = REFERENCE_GENOME
#     shell:
#         '''
#         scripts/mosdepth_to_x_coverage.py \
#         -f {params.genome} \
#         -m {input.bed} \
#         > {output}
#         '''