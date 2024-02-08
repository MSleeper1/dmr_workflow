
rule multiqc_dir:
    input:
        expand("samtools_stats/{sample}.txt", sample=["a", "b"]),
    output:
        "qc/multiqc.html",
        directory("qc/multiqc_data"),
    params:
        extra="--data-dir",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.3.6/bio/multiqc"