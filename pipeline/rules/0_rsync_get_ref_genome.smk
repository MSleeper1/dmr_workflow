
rule rsync_get_ref_genome:
    output:
        zipped_fasta = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),
        unzipped_fasta = expand("{ref_dir}/{fasta}.fa", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])
    
    log:
        expand("../pre-processing/logs/rule-logs/rsync_get_ref_genome/rsync_get_ref_genome--{fasta}.log", fasta = config["ref"]["fasta"])

    conda:
        "../env/rsync.yaml"
    
    params:
        gold_path = config["ref"]["goldenPath"], 
        out_dir = config["ref"]["dir"],
        rsync_args = config["prep_args"]["rsync_get_ref"],
        unzipped_fasta = expand("{ref_dir}/{fasta}.fa", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])

    shell:
        """
        mkdir -p {params.out_dir}
        rsync {params.rsync_args} rsync:{params.gold_path} {params.out_dir} > {log} 2>&1
        gunzip -c {output.zipped_fasta} > {output.unzipped_fasta} 
        """