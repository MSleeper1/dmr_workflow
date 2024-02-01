

rule wgbstools_init_ref_genome:
    input:
        ref = expand("{ref_dir}/{fasta}.fa.gz", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"])

    output:
        directory(expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome=config["ref"]["wgbstools_ref_name"]))
    
    log:
        stdout = expand("../logs/wgbstools_init_ref_genome/wgbstools_init_ref_genome--{fasta}.out", fasta=config["ref"]["fasta"]),
        stderr = expand("../logs/wgbstools_init_ref_genome/wgbstools_init_ref_genome--{fasta}.err", fasta=config["ref"]["fasta"])

    conda:
        "../env/wgbstools.yaml"

    params:
        genome=config["ref"]["wgbstools_ref_name"],
        unzipped_ref = expand("{ref_dir}/{fasta}.fa", ref_dir=config["ref"]["dir"], fasta=config["ref"]["fasta"]),

    shell:
        """
        wgbstools init_genome {params.genome} --fasta_path {params.unzipped_ref} > {log.stdout} 2> {log.stderr}
        """
