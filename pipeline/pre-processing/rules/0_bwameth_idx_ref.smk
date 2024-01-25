

rule index_bwameth_ref_genome:
    input:
        ref = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),
    
    output:
        expand("{ref_dir}/{fasta}.fa.gz.bwameth.{suf}", suf=["c2t","c2t.bwt","c2t.pac","c2t.ann","c2t.amb","c2t.sa"], ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),

    log:
        stdout="logs/0_index_bwameth_ref.stdout",
        stderr="logs/0_index_bwameth_ref.stderr",

    conda:
        "../env/bwameth.yaml",

    shell:
        "bwameth.py index {input.ref} > {log.stdout} 2> {log.stderr}"