### Genome reference file download rules ###

# rsync_get_ref_genome: download the reference genome from the ref goldenPath specified in the config file
# change ref.goldenPath and ref.fasta in the config file to the desired goldenPath and fasta file name expected
rule rsync_get_ref_genome:
    output:
        zipped_fasta = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),
        unzipped_fasta = expand("{ref_dir}/{fasta}.fa", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])
    
    log:
        expand("../pre-processing/logs/rule-logs/00_genome_prep/00_rsync_get_ref_genome--{fasta}.log", fasta = config["ref"]["fasta"])

    conda:
        "../environment_files/rsync.yaml"
    
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

# rsync_get_ref_gtf: download the reference gtf from the gtf goldenPath specified in the config file
# change ref.gtf_goldenPath and ref.gtf in the config file to the desired goldenPath and gtf file name expected
rule rsync_get_ref_gtf:
    output:
        zipped_gtf = expand("{ref_dir}/{gtf}.gtf.gz", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"]),
        unzipped_gtf = expand("{ref_dir}/{gtf}.gtf", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])
    
    log:
        expand("../pre-processing/logs/rule-logs/rsync_get_ref_gtf/rsync_get_ref_gtf--{gtf}.log", gtf = config["ref"]["gtf"])

    conda:
        "../environment_files/rsync.yaml"
    
    params:
        gold_path = config["ref"]["gtf_goldenPath"], 
        out_dir = config["ref"]["dir"],
        rsync_args = config["prep_args"]["rsync_get_ref"]

    shell:
        """
        mkdir -p {params.out_dir}
        rsync {params.rsync_args} rsync:{params.gold_path} {params.out_dir} > {log} 2>&1
        gunzip -c {output.zipped_gtf} > {output.unzipped_gtf}
        """