### Genome preparation for bismark ###

# bismak_index_ref_genome: index the reference genome for bismark
# change ref.dir and ref.fasta in config.yaml to use a different reference genome
rule bismark_index_ref_genome:
    input:
        fasta_path = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])

    output:
        directory(expand("{bismark_idx_dir}", bismark_idx_dir=config["ref"]["bismark_idx_dir"]))
    
    log:
        stdout = expand("../pre-processing/logs/rule-logs/00_genome_prep/00_bismark_index_ref_genome--{fasta}.out", fasta = config["ref"]["fasta"]),
        stderr = expand("../pre-processing/logs/rule-logs/00_genome_prep/00_bismark_index_ref_genome--{fasta}.err", fasta = config["ref"]["fasta"])
    
    conda:
        "../environment_files/bismark.yaml"
    
    shadow:
        "shallow"

    params:
        bismk_args = config["prep_args"]["bismark_genome_prep"],
        bismk_dir = expand("{ref_dir}/bismark", ref_dir = config["ref"]["dir"]),
        bismk_fasta_path = expand("{ref_dir}/bismark/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])

    shell:
        """
        mkdir -p {params.bismk_dir} 
        cp {input.fasta_path} {params.bismk_fasta_path} 
        bismark_genome_preparation {params.bismk_args} {params.bismk_dir} > {log.stdout} 2> {log.stderr}
        """


