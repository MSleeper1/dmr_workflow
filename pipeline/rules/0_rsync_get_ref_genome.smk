
rule rsync_get_ref_genome:
    output:
        zipped_fasta = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),
        unzipped_fasta = expand("{ref_dir}/{fasta}.fa", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])
    
    log:
        stdout = expand("../logs/rsync_get_ref_genome/rsync_get_ref_genome--{fasta}.out", fasta = config["ref"]["fasta"]),
        stderr = expand("../logs/rsync_get_ref_genome/rsync_get_ref_genome--{fasta}.err", fasta = config["ref"]["fasta"])

    conda:
        "../env/rsync.yaml"
    
    params:
        gold_path = config["ref"]["goldenPath"], 
        out_dir = config["ref"]["dir"],
        rsync_args = config["prep_args"]["rsync_get_ref"]

    shell:
        """
        mkdir -p {params.out_dir} > {log.stdout} 2> {log.stderr}
        rsync {params.rsync_args} rsync:{params.gold_path} {params.out_dir} > {log.stdout} 2> {log.stderr}
        gunzip -c {output.zipped_fasta} > {output.unzipped_fasta} > {log.stdout} 2> {log.stderr}
        """
