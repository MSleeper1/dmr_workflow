### Genome preparation for fastq_screen ###

# fastq_screen_genome_prep: rule to download the fastq_screen genomes that have been already been prepared with bismark genome preparation
# This rule will download the fastq_screen genomes that have been already been prepared with bismark genome preparation
# Fastq_screen will create the directory and download the genomes to the directory and create a fastq_screen.conf file
# The fastq_screen.conf file will be used to specify the genomes to be used in the fastq_screen rule
# default genomes included by fastqc are human, mouse, rat, drosophilia, worm, yeast, arabdopsis, ecoli, and PhiX
# the default genomes can be changed by editing the fastq_screen.conf file
# An alternative method is to download the genomes of interest from ncbi and prepare them with bismark_genome_preparation

rule fastq_screen_genome_prep:
    output:
        directory(expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/", genomes_dir=config["genomes_dir"])),
        expand("{genomes_dir}/FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf", genomes_dir=config["genomes_dir"])
    
    log:
        "../pre-processing/logs/rule-logs/00_genome_prep/00_fastq_screen_genome_prep.err"

    conda:
        "../environment_files/fastq-screen.yaml"

    # shadow: 
    #     "shallow"

    params:
        genomes_dir=expand("{genomes_dir}", genomes_dir=config["genomes_dir"])
    
    shell:
        """
        mkdir -p {params.genomes_dir}
        cd {params.genomes_dir}
        fastq_screen --bisulfite --get_genomes >2 {log}
        """