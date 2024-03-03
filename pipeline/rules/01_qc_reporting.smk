### Snakemake rules for quality control reporting of raw sequence files (01) ###
# fastqc and fastq_screen are used to generate quality control reports for raw sequence files

### FASTQC ###
# FastQC is a quality control tool for high throughput sequence data. It reads in sequence data in a variety of formats and will create an HTML report to view the results.
# input: raw sequence files (fastq)
# output: fastqc report (html, zip)

# Rule to run fastqc on single-end sequence files
rule fastqc_se:
    input:
        expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir = config["data"]["dir"])
    output:
        expand("{rep_dir}/01_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])
    log:
        "../pre-processing/logs/rule-logs/01_fastqc_se/{ref}/01_fastqc_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
    conda:
        "../env/fastqc.yaml"
    params:
        output_dir = expand("{rep_dir}/01_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir = config["reports_dir"])
    wildcard_constraints:
        layout = "se"
    shell:
        """
        mkdir -p {params.output_dir}
        fastqc -o {params.output_dir} {input} > {log} 2>&1
        """

# Rule to run fastqc on paired-end sequence files
rule fastqc_pe:
	input: 
		r1 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", data_dir = config["data"]["dir"]),
		r2 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", data_dir = config["data"]["dir"])
	output:
		r1 = expand("{rep_dir}/01_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"]),
		r2 = expand("{rep_dir}/01_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"])
	log:
        "../pre-processing/logs/rule-logs/01_fastqc_pe/{ref}/01_fastqc_pe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.log"
	conda:
		"../env/fastqc.yaml"
	params:
		output_dir = expand("{rep_dir}/01_fastqc/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir = config["reports_dir"])
	wildcard_constraints:
		layout = "pe"
	shell: 
		"""
		mkdir -p {params.output_dir}
		fastqc -o {params.output_dir} {input.r1} {input.r2} > {log} 2>&1
		"""

### FASTQ SCREEN ###
# FastQ Screen is a tool to screen a set of sequences in FastQ format against a set of sequence databases so you can see if the composition of the library matches with what you expect.
# input: raw sequence files (fastq)
# output: fastq_screen report (txt, html)

# Rule to run fastq_screen on single-end sequence files
rule fastq_screen_se:
    input:
        conf = expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", genomes_dir=config["genomes_dir"]),
        fq_file = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir=config["data"]["dir"])
    output:
        expand("{rep_dir}/01_fastq_screen/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_screen.{suf}", rep_dir=config["reports_dir"], suf=["txt", "html"])
    log:
        "../pre-processing/logs/rule-logs/01_fastq_screen_se/{ref}/01_fastq_screen_se-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"
    conda:
        "../env/fastq-screen.yaml"
    params:
        out_dir = expand("{rep_dir}/01_fastq_screen/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"])
    threads: 6
    wildcard_constraints:
        layout = "se"
    shell:
        """
        mkdir -p {params.out_dir}
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.fq_file} >2 {log}
        """

# Rule to run fastq_screen on paired-end sequence files
rule fastq_screen_pe:
    input:
        conf = expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", genomes_dir=config["genomes_dir"]),
        r1 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1.fastq", data_dir = config["data"]["dir"]),
        r2 = expand("{data_dir}/01_raw_sequence_files/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2.fastq", data_dir = config["data"]["dir"])
    output:
        expand("{rep_dir}/01_fastq_screen/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_1_screen.{suf}", rep_dir=config["reports_dir"], suf=["txt", "html"]),
        expand("{rep_dir}/01_fastq_screen/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}/{{accession}}_2_screen.{suf}", rep_dir=config["reports_dir"], suf=["txt", "html"])
    log:
        "../pre-processing/logs/rule-logs/01_fastq_screen_pe/{ref}/01_fastq_screen_pe-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"
    conda:
        "../env/fastq-screen.yaml"
    threads: 6
    wildcard_constraints:
        layout = "pe"
    shell:
        """
        mkdir -p {params.out_dir}
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.r1} >2 {log}
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.r2} >2 {log}
        """