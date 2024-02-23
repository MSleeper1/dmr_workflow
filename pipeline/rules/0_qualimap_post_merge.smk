### Qualimap of merged bam files ###
# Qualimaps flag --outdir is a relative path. To avoid issues with this, the ouput will be moved to the reports directory after qualimap is finished.

# place in snakefile
#   ### qualimap output report directories ###
#   expand("{rep_dir}/qualimap/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output

rule qualimap_post_merge:
    input:
        expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"])

    output:
        directory(expand("{rep_dir}/qualimap/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged", rep_dir=config["reports_dir"]))
        
    log:
        "../pre-processing/logs/rule-logs/qualimap_post_merge/{ref}/qualimap_post_merge-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    conda:
        "../env/qualimap.yaml"

    params:
        temp_out = "{ref}-{patient_id}-{group}-{srx_id}-{layout}_merged"

    shell:
        """
        qualimap bamqc -bam {input} -c -sd -os -gd hg38 -gff /home/msleeper/scratch/genomes/hg38/analysisSet/hg38.knownGene.gtf --outdir {params.temp_out} > {log} 2>&1
        mv -f -v {params.temp_out} {output} >> {log} 2>&1
        rm -rf -v {params.temp_out} >> {log} 2>&1
        """