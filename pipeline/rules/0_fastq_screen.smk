

rule fastq_screen:
    input:
        conf = expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", genomes_dir=config["genomes_dir"]),
        fq_file = expand("{data_dir}/raw_sequence_files/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}.fastq", data_dir=config["data"]["dir"])
        
    output:
        expand("{rep_dir}/fastq_screen/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}/{{accession}}_screen.{suf}", rep_dir=config["reports_dir"], suf=["txt"])

    log:
        "../pre-processing/logs/rule-logs/fastq_screen/{ref}/fastq_screen-{ref}-{patient_id}-{group}-{srx_id}-{layout}-{accession}.err"

    conda:
        "../env/fastq-screen.yaml"

    params:
        out_dir = expand("{rep_dir}/fastq_screen/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}-{{layout}}", rep_dir=config["reports_dir"])

    threads: 6

    shell:
        """
        mkdir -p {params.out_dir}
        fastq_screen --bisulfite --aligner bowtie2 --conf {input.conf} --threads {threads} --outdir {params.out_dir} --force --quiet {input.fq_file} >2 {log}
        """
