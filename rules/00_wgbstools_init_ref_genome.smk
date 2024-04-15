### Genome prepatation for wgbstools use ###

# wgbstools_init_ref_genome: Initialize reference genome for wgbstools
# change ref.fasta and ref.wgbstools_ref_name in config.yaml to the name of the desired reference genome fasta file and name for wgbstools to use for the reference genome
rule wgbstools_init_ref_genome:
    input:
        ref = expand("{ref_dir}/{fasta}.fa.gz", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"])

    output:
        directory(expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"]))
    
    log:
        expand("../pre-processing/rule-logs/00_genome_prep/00_wgbstools_init_ref_genome--{fasta}.log", fasta=config["ref"]["fasta"])

    conda:
        "../environment_files/wgbstools.yaml"

    shadow:
        "shallow"

    params:
        genome=config["ref"]["wgbstools_ref_name"],
        unzipped_ref = expand("{ref_dir}/{fasta}.fa", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),

    shell:
        """
        wgbstools init_genome {params.genome} --fasta_path {params.unzipped_ref} > {log} 2>&1
        """
