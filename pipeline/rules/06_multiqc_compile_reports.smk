
rule multiqc_compile_reports:
    input:
        expand("{rep_dir}/01_fastqc/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{rep_dir}/01_fastqc/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{rep_dir}/01_fastq_screen/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{rep_dir}/01_fastq_screen/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{rep_dir}/02_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
        expand("{rep_dir}/02_fastqc_post_trim/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
        expand("{rep_dir}/02_fastqc_post_trim/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
        expand("{rep_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim
        expand("{rep_dir}/03_bwameth/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bwameth_mapping se report from stderr
        expand("{rep_dir}/03_bwameth/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bwameth_mapping pe report from stderr
        expand("{rep_dir}/03_bismark_bwt2/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bismark_mapping se report
        expand("{rep_dir}/03_bismark_bwt2/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bismark_mapping pe report
        expand("{rep_dir}/04_fastqc_post_dedup/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc se output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output
        expand("{rep_dir}/05_fastqc_post_merge/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc post-merge output
        expand("{rep_dir}/05_feature_counts/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample=sample_info.itertuples()), # featureCounts output
        expand("{rep_dir}/05_qualimap/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output
        expand("{rep_dir}/05_samtools_post_merge/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output
    output:
        expand("{rep_dir}/prep_multiqc_data/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/prep_multiqc_data", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc.log"
    conda:
        "../env/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"