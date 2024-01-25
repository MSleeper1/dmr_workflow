
rule rsync_get_ref:
    output:
        expand("{outdir}/{fasta}.fa.gz", outdir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),
    params:
        # set params.path equal to the goldenPath of the reference genome in the config file
        path = config["ref"]["goldenPath"],
        
        # set params.outputdir equal to the output directory in the config file
        outputdir = config["ref"]["dir"],

    log:
        stdout = "logs/rsync_get_ref.stdout",
        stderr = "logs/rsync_get_ref.stderr",

    conda:
        "../env/rsync.yaml",

    shell:
        "rsync -avzP --progress --ignore-existing rsync:{params.path} {params.outputdir} > {log.stdout} 2> {log.stderr}"
