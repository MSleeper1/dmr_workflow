
#outputs not added to snakefile yet
#add inport of rule and expected outputs to snakefile for piped version to test
# outputs need to be created for se and pe versions of the piped rule

rule bwameth_mapping_se:
	input:
        read_se = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa", bwa_idx_dir=config["ref"]["bwa_idx_dir"], ref=config["ref"]["fasta"])

	output:
		sam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.sam", data_dir=config["data"]["dir"]),
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.bam", data_dir=config["data"]["dir"])   

	threads: 4

	log:
		"logs/mapping/bwameth/{tool}/{trim_galoreparameterSet}/{species}/{sample}.log"

	conda:
		"../envs/bwameth.yaml"

    params:
        accession = "{accession}"

	shell:
        """
        echo "aligning {input.read_se} to {input.index}"
        bwameth.py --threads {threads} --reference {input.index} {input.read_se} > {output.sam}
        echo "converting {output.sam} to {output.bam}"
        samtools view -S -b {output.sam} > {output.bam}
        echo "alignment and conversion complete for {accession}"
        """

rule bwameth_mapping_se_pipe:
	input:
        read_se = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa", bwa_idx_dir=config["ref"]["bwa_idx_dir"], ref=config["ref"]["fasta"])

	output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.bam", data_dir=config["data"]["dir"]) 

	threads: 4

	log:
		"logs/mapping/bwameth/{tool}/{trim_galoreparameterSet}/{species}/{sample}.log"

	conda:
		"../envs/bwameth.yaml"

    params:
        accession = "{accession}

    shell:
        """
        echo "aligning {input.read_se} to {input.index} and converting sam output to {output.bam}"
        bwameth.py --threads {threads} --reference {input.index} {input.read_se} | samtools view -S -b > {output.bam}
        echo "bwameth single end alignment and sam to bam conversion complete for {accession}"
        """

rule bwameth_mapping_pe_pipe:
	input:
        read_1 = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_1.fq", data_dir=config["data"]["dir"]),
        read_2 = expand("{data_dir}/trimmed/trim_galore/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}_2.fq", data_dir=config["data"]["dir"]),
        index = expand("{bwa_idx_dir}/{fasta}.fa", bwa_idx_dir=config["ref"]["bwa_idx_dir"], ref=config["ref"]["fasta"])

	output:
        bam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/{{ref}}/{{patient_id}}/{{group}}-{{srx_id}}/{{accession}}.bam", data_dir=config["data"]["dir"]) 

	threads: 4

	log:
		"logs/mapping/bwameth/{tool}/{trim_galoreparameterSet}/{species}/{sample}.log"

	conda:
		"../envs/bwameth.yaml"

    params:
        accession = "{accession}

    shell:
        """
        echo "aligning {input.read_1} and {input.read_2} to {input.index} and converting sam output to {output.bam}"
        bwameth.py --threads {threads} --reference {input.index} {input.read_1} {input.read_2} | samtools view -S -b > {output.bam}
        echo "bwameth paired end alignment and sam to bam conversion complete for {accession}"
        """