# snakemake rule to index reference genome for bismark bwt2 alignment

rule index_bismark_bwt2_ref_genome:
    input:
        ref_dir = config["ref"]["dir"],

    output:
        directory(expand("{ref_dir}/Bisulfite_Genome", ref_dir=config["ref"]["dir"])),
    
    log:
        stdout="logs/0_index_bismark_ref.stdout",
        stderr="logs/0_index_bismark_ref.stderr",

    conda:
        "../env/bismark.yaml",

    params:
        bismk = config["prep_args"]["bismark_genome_prep"]

    shell:
        "bismark_genome_preparation {params.bismk} {input.ref_dir} > {log.stdout} 2> {log.stderr}"

    

