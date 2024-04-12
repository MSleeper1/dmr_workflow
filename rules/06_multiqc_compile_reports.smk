
### Multiqc rules for compiling reports from pre-processing steps

## compile all reports from pre-processing steps
# compile reports for bwameth mapped files at each step
rule multiqc_compile_reports_bwa:
    input:
        expand("{rep_dir}/01_fastqc/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{rep_dir}/01_fastqc/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{rep_dir}/01_fastq_screen/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{rep_dir}/01_fastq_screen/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        # 02 trimming reports from trim_galore and fastqc for pe and se layouts
        expand("{rep_dir}/02_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
        expand("{rep_dir}/02_fastqc_post_trim/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
        expand("{rep_dir}/02_fastqc_post_trim/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
        expand("{rep_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim
        # 03 mapping reports from bwameth 
        expand("{rep_dir}/03_bwameth/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], sample=sample_info_se.itertuples()), # bwameth_mapping reports from stderr
        # 04 dedup reports from sambamba
        expand("{rep_dir}/04_sambamba/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # sambamba dedup report log for bwa mapped reads
        # 04 QC post-dedup 
        expand("{rep_dir}/04_fastqc_post_dedup_bwa/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bwa output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bwa/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bwameth mapped
        # 05 QC post merge 
        expand("{rep_dir}/05_fastqc_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc post-merge output
        expand("{rep_dir}/05_feature_counts_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample=sample_info.itertuples()), # featureCounts output
        expand("{rep_dir}/05_qualimap_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output for bwameth mapped reads
        expand("{rep_dir}/05_samtools_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output       
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.global.dist.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.per-base.bed.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # produced unless --no-per-base specified in mosdepth rule
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.summary.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
    output:
        expand("{rep_dir}/06_multiqc_data/bwa/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/bwa", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_bwa.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# compile reports for bismark mapped files at each step
rule multiqc_compile_reports_bis:
    input:
        expand("{rep_dir}/01_fastqc/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{rep_dir}/01_fastqc/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{rep_dir}/01_fastq_screen/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{rep_dir}/01_fastq_screen/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        # 02 trimming reports from trim_galore and fastqc for pe and se layouts
        expand("{rep_dir}/02_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
        expand("{rep_dir}/02_fastqc_post_trim/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
        expand("{rep_dir}/02_fastqc_post_trim/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
        expand("{rep_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim
        # 03 mapping reports from bismark for pe and se layouts
        expand("{rep_dir}/03_bismark_bwt2/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bismark_mapping se report
        expand("{rep_dir}/03_bismark_bwt2/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bismark_mapping pe report
        # 04 dedup reports from bismark
        expand("{rep_dir}/04_bismark_deduplication/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark_deduplication_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # bismark_deduplication output report for bismark mapped reads
        # 04 QC post-dedup
        expand("{rep_dir}/04_fastqc_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_deduped_bismark_bt2_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bis output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_deduped_bismark_bt2.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bismark mapped
        # 05 QC post merge
        expand("{rep_dir}/05_fastqc_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc post-merge output
        expand("{rep_dir}/05_feature_counts_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample=sample_info.itertuples()), # featureCounts output
        expand("{rep_dir}/05_qualimap_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output for bismark mapped reads
        expand("{rep_dir}/05_samtools_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.global.dist.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.per-base.bed.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # produced unless --no-per-base specified in mosdepth rule
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.summary.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        # 06 methylation extraction reports
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.M-bias.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_splitting_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bismark.cov.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bedGraph.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"])
    output:
        expand("{rep_dir}/06_multiqc_data/bis/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/bis", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_bis.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

## compile reports by pre-processing step
# 01 raw read reports
rule multiqc_compile_reports_01:
    input:
        expand("{rep_dir}/01_fastqc/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{rep_dir}/01_fastqc/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{rep_dir}/01_fastq_screen/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{rep_dir}/01_fastq_screen/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
    output:
        expand("{rep_dir}/06_multiqc_data/01_raw/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/01_raw", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_01.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 02 trimmed read reports
rule multiqc_compile_reports_02:
    input:
        expand("{rep_dir}/02_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
        expand("{rep_dir}/02_fastqc_post_trim/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
        expand("{rep_dir}/02_fastqc_post_trim/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
        expand("{rep_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim
    output:
        expand("{rep_dir}/06_multiqc_data/02_trimmed/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/02_trimmed", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_02.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 03 bwameth mapping reports
rule multiqc_compile_reports_03_bwa:
    input:
        expand("{rep_dir}/03_bwameth/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], sample=sample_info_se.itertuples()), # bwameth_mapping reports from stderr
    output:
        expand("{rep_dir}/06_multiqc_data/03_bwameth_mapped/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/03_bwameth_mapped", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_03_bwa.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 03 bismark mapping reports
rule multiqc_compile_reports_03_bis:
    input:
        expand("{rep_dir}/03_bismark_bwt2/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bismark_mapping se report
        expand("{rep_dir}/03_bismark_bwt2/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bismark_mapping pe report
    output:
        expand("{rep_dir}/06_multiqc_data/03_bismark_mapped/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/03_bismark_mapped", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_03_bis.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 04 bwameth mapped: post deduplication
rule multiqc_compile_reports_04_bwa:
    input:
        expand("{rep_dir}/04_sambamba/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # sambamba dedup report log for bwa mapped reads
        expand("{rep_dir}/04_fastqc_post_dedup_bwa/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bwa output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bwa/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bwameth mapped
    output:
        expand("{rep_dir}/06_multiqc_data/04_sambamba_deduped_bwa/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/04_sambamba_deduped_bwa", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_04_bwa.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 04 bismark mapped: post deduplication
rule multiqc_compile_reports_04_bis:
    input:
        expand("{rep_dir}/04_bismark_deduplication/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark_deduplication_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # bismark_deduplication output report for bismark mapped reads
        expand("{rep_dir}/04_fastqc_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_deduped_bismark_bt2_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bis output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_deduped_bismark_bt2.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bismark mapped
    output:
        expand("{rep_dir}/06_multiqc_data/04_bismark_deduped_bis/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/04_bismark_deduped_bis", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_04_bis.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 05 bwameth mapped: post-merge
rule multiqc_compile_reports_05_bwa:
    input:
        expand("{rep_dir}/05_fastqc_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc post-merge output
        expand("{rep_dir}/05_feature_counts_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample=sample_info.itertuples()), # featureCounts output
        expand("{rep_dir}/05_qualimap_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output for bwameth mapped reads
        expand("{rep_dir}/05_samtools_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output       
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.global.dist.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.per-base.bed.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # produced unless --no-per-base specified in mosdepth rule
        expand("{rep_dir}/05_mosdepth_post_merge_bwa/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.summary.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
    output:
        expand("{rep_dir}/06_multiqc_data/05_post_merge_bwa/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/05_post_merge_bwa", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_05_bwa.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 05 bismark mapped: post-merge
rule multiqc_compile_reports_05_bis:
    input:
        expand("{rep_dir}/05_fastqc_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc post-merge output
        expand("{rep_dir}/05_feature_counts_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.featureCounts{suf}", rep_dir=config["reports_dir"], suf=["", ".summary", ".jcounts"], sample=sample_info.itertuples()), # featureCounts output
        expand("{rep_dir}/05_qualimap_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # qualimap output for bismark mapped reads
        expand("{rep_dir}/05_samtools_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.global.dist.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.per-base.bed.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # produced unless --no-per-base specified in mosdepth rule
        expand("{rep_dir}/05_mosdepth_post_merge_bis/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_merged.mosdepth.summary.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"])
    output:
        expand("{rep_dir}/06_multiqc_data/05_post_merge_bis/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/05_post_merge_bis", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_05_bis.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"

# 06 bismark mapped: methylation extraction reports
rule multiqc_compile_reports_06_bis:
    input:
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.M-bias.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_splitting_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bismark.cov.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bedGraph.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"])
    output:
        expand("{rep_dir}/06_multiqc_data/06_methyl_extract_bis/multiqc.html", rep_dir=config["reports_dir"]),
        directory(expand("{rep_dir}/06_multiqc_data/06_methyl_extract_bis", rep_dir=config["reports_dir"]))
    log:
        "../pre-processing/logs/rule-logs/06_prep_multiqc_06_methyl_bis.log"
    conda:
        "../environment_files/wgbstools.yaml"
    params:
        extra="--verbose --fullnames"  # Optional: extra parameters for multiqc.
    wrapper:
        "v3.4.1/bio/multiqc"
