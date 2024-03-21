### Genome preparation for bwameth ###

# bwameth_index_ref_genome: index the reference genome for bwameth
# change ref.dir and ref.fasta in config.yaml to use a different reference genome
rule bwameth_index_ref_genome:
    input:
        fasta_path = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])
    
    output:
        expand("{bwa_idx_dir}/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t","c2t.bwt","c2t.pac","c2t.ann","c2t.amb","c2t.sa"], bwa_idx_dir = config["ref"]["bwa_idx_dir"], fasta = config["ref"]["fasta"]),
        expand("{bwa_idx_dir}/{fasta}.fa.gz", bwa_idx_dir = config["ref"]["bwa_idx_dir"], fasta = config["ref"]["fasta"])

    log:
        expand("../pre-processing/logs/rule-logs/bwameth_index_ref_genome/bwameth_index_ref_genome--{fasta}.log", fasta = config["ref"]["fasta"])

    conda:
        "../environment_files/bwameth.yaml"

    params:
        bwa_dir = config["ref"]["bwa_idx_dir"],
        fasta = config["ref"]["fasta"]

    shell:
        """
        mkdir -p {params.bwa_dir} 
        cp {input.fasta_path} {params.bwa_dir}/{params.fasta}.fa.gz
        bwameth.py index {params.bwa_dir}/{params.fasta}.fa.gz > {log} 2>&1
        """