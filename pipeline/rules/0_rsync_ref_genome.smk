rule rsync_ref_genome:
    output:
        expand("{outdir}{fasta}", outdir = config["ref"]["outputdir"], fasta = config["ref"]["fasta"]),

    params:
        # set params.path equal to the goldenPath of the reference genome in the config file
        path = config["ref"]["goldenPath"],

        # set params.outputdir equal to the output directory in the config file
        outputdir = config["ref"]["outputdir"],
        
    log: 
        "logs/rsync_ref_genome.log",
    
    conda: 
        "../env/rsync.yaml"
    
    shell: 
        "rsync -avzP --progress --ignore-existing rsync:{params.path} {params.outputdir}"