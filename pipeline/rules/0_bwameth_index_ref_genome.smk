

rule bwameth_index_ref_genome:
    input:
        fasta_path = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])
    
    output:
        expand("{bwa_idx_dir}/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t","c2t.bwt","c2t.pac","c2t.ann","c2t.amb","c2t.sa"], bwa_idx_dir = config["ref"]["bwa_idx_dir"], fasta = config["ref"]["fasta"])

    log:
        stdout = expand("../logs/bwameth_index_ref_genome/bwameth_index_ref_genome--{fasta}.out", fasta = config["ref"]["fasta"]),
        stderr = expand("../logs/bwameth_index_ref_genome/bwameth_index_ref_genome--{fasta}.err", fasta = config["ref"]["fasta"])

    conda:
        "../env/bwameth.yaml"

    params:
        bwa_dir = config["ref"]["bwa_idx_dir"],
        fasta = config["ref"]["fasta"]

    shell:
        """
        mkdir -p {params.bwa_dir} > {log.stdout} 2> {log.stderr}
        cp {input.fasta_path} {params.bwa_dir}/{params.fasta}.fa.gz > {log.stdout} 2> {log.stderr}
        bwameth.py index {params.bwa_dir}/{params.fasta}.fa.gz > {log.stdout} 2> {log.stderr}
        """