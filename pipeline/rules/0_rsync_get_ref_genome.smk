
rule rsync_get_ref:
    output:
        expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),

    params:
        # set params.path equal to the goldenPath of the reference genome in the config file
        path = config["ref"]["goldenPath"],
        
        # set params.outputdir equal to the reference directory in the config file
        outputdir = config["ref"]["dir"],

        rsync = config["prep_args"]["rsync_get_ref"],

    log:
        stdout = "logs/rsync_get_ref.stdout",
        stderr = "logs/rsync_get_ref.stderr",

    conda:
        "../env/rsync.yaml",

    shell:
        "rsync {params.rsync} rsync:{params.path} {params.outputdir} > {log.stdout} 2> {log.stderr}"
