# snakemake rule to index reference genome for bismark bwt2 alignment

rule bismark_index_ref_genome:
    input:
        fasta_path = expand("{ref_dir}/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"])

    output:
        directory(expand("{bismark_idx_dir}", bismark_idx_dir=config["ref"]["bismark_idx_dir"]))
    
    log:
        stdout = expand("../logs/bismark_index_ref_genome/bismark_index_ref_genome--{fasta}.out", fasta = config["ref"]["fasta"]),
        stderr = expand("../logs/bismark_index_ref_genome/bismark_index_ref_genome--{fasta}.err", fasta = config["ref"]["fasta"])
    
    conda:
        "../env/bismark.yaml"

    params:
        bismk_args = config["prep_args"]["bismark_genome_prep"],
        bismk_dir = expand("{ref_dir}/bismark", ref_dir = config["ref"]["dir"]),
        bismk_fasta_path = expand("{ref_dir}/bismark/{fasta}.fa.gz", ref_dir = config["ref"]["dir"], fasta = config["ref"]["fasta"]),


    shell:
        """
        mkdir -p {params.bismk_dir} 
        cp {input.fasta_path} {params.bismk_fasta_path} 
        bismark_genome_preparation {params.bismk_args} {params.bismk_dir} > {log.stdout} 2> {log.stderr}
        """

    

