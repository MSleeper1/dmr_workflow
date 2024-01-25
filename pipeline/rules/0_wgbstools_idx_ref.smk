

rule index_wgbstools_ref_genome:
    input:
        ref = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta=config["ref"]["fasta"]),

    output:
        directory(expand("{wgbstools_ref_dir}/{genome}", wgbstools_ref_dir=config["ref"]["wgbstools_idx_dir"], genome = config["ref"]["wgbstools_ref_name"])),
    
    log:
        stdout="logs/index_wgbstools_ref.stdout",
        stderr="logs/index_wgbstools_ref.stderr",

    conda:
        "../env/wgbstools.yaml",

    params:
        genome = config["ref"]["wgbstools_ref_name"],

    shell:
        "wgbstools init_genome {params.genome} --fasta_path {input.ref} > {log.stdout} 2> {log.stderr}"