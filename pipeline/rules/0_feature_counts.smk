rule feature_counts:
    input:
        sam = expand("{data_dir}/trimmed/trim_galore/aligned/bwameth/deduped/sambamba/merged/{{ref}}-{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.bam", data_dir=config["data"]["dir"]),
        annotation = expand("{ref_dir}/{gtf}", ref_dir = config["ref"]["dir"], gtf = config["ref"]["gtf"])
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    
    output:
        expand("{rep_dir}/feature_counts/{{ref}}/{{patient_id}}-{{group}}-{{srx_id}}-{{layout}}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"])

    log:
        "../pre-processing/logs/rule-logs/feature_counts/{ref}/feature_counts-{ref}-{patient_id}-{group}-{srx_id}-{layout}.log"

    threads:
        3

    conda:
        "../env/feature_counts.yaml"

    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -f"

    wrapper:
        "0.72.0/bio/subread/featurecounts"