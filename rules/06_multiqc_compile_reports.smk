
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
        expand("{rep_dir}/04_sambamba_dedup/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # sambamba dedup report log for bwa mapped reads
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
        file = expand("{rep_dir}/summary/01-06_bwa_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/01-06_bwa_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/01-06_bwa_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "01-06_bwa_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for all pre-processing steps for bwameth mapped reads. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "01-06 All QC Reports for bwameth mapped files", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for all qc reports in bwameth mapped workflow: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

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
        expand("{rep_dir}/04_bismark_deduplication/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}__bismark.deduplication_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # bismark_deduplication output report for bismark mapped reads
        # 04 QC post-dedup
        expand("{rep_dir}/04_fastqc_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bis output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bismark mapped
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
        file = expand("{rep_dir}/summary/01-06_bis_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/01-06_bis_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/01-06_bis_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "01-06_bis_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for all pre-processing steps for bismark mapped reads. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "01-06 All QC Reports for bismark mapped files", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for all qc reports in bwameth mapped workflow: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

## compile reports by pre-processing step
# 01 raw read reports
rule multiqc_compile_reports_01:
    input:
        expand("{rep_dir}/01_fastqc/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # se fastqc se output
        expand("{rep_dir}/01_fastqc/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html", "zip"]), # pe fastqc pe R1 and R2 output
        expand("{rep_dir}/01_fastq_screen/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_screen.{suf}", se=sample_info_se.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"]), # se fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
        expand("{rep_dir}/01_fastq_screen/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_screen.{suf}", pe=sample_info_pe.itertuples(), suf=["txt", "html"], rep_dir=config["reports_dir"], read=["1", "2"]), # pe fastq_screen output (other outputs: "png", "html", "bisulfite_orientation.png")
    output:
        file = expand("{rep_dir}/summary/01_raw_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/01_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/01_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "01_raw_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for raw reads produced by fastqc and fastq_screen. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "01 Raw Sequence QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """


# 02 trimmed read reports
rule multiqc_compile_reports_02:
    input:
        expand("{rep_dir}/02_trim_galore/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # moved trim_galore se trimming report
        expand("{rep_dir}/02_fastqc_post_trim/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], se=sample_info_se.itertuples(), suf=["html","zip"]), # fastqc report post-trim
        expand("{rep_dir}/02_fastqc_post_trim/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}_trimmed_fastqc.{suf}", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"], suf=["html","zip"]), # trimmed fastq outputs
        expand("{rep_dir}/02_trim_galore/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_{read}.fastq_trimming_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples(), read=["1", "2"]), # fastqc reports post-trim
    output:
        file = expand("{rep_dir}/summary/02_trimmed_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/02_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/02_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "02_trimmed_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for trimmed reads produced by trim_galore and fastqc. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "02 Trimmed read QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

# 03 bwameth mapping reports
rule multiqc_compile_reports_03_bwa:
    input:
        expand("{rep_dir}/03_bwameth/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_bwameth_report.txt", rep_dir=config["reports_dir"], sample=sample_info_se.itertuples()), # bwameth_mapping reports from stderr
    output:
        file = expand("{rep_dir}/summary/03_bwameth_mapping_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/03_bwa_map_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/03_bwa_map_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "03_bwameth_mapping_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for reads mapped by bwameth. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "03 Bwameth Mapping QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

# 03 bismark mapping reports
rule multiqc_compile_reports_03_bis:
    input:
        expand("{rep_dir}/03_bismark_bwt2/{se.ref}/{se.patient_id}-{se.group}-{se.srx_id}-{se.layout}/{se.accession}_trimmed_bismark_bt2_SE_report.txt", rep_dir=config["reports_dir"], se=sample_info_se.itertuples()), # bismark_mapping se report
        expand("{rep_dir}/03_bismark_bwt2/{pe.ref}/{pe.patient_id}-{pe.group}-{pe.srx_id}-{pe.layout}/{pe.accession}_trimmed_bismark_bt2_PE_report.txt", rep_dir=config["reports_dir"], pe=sample_info_pe.itertuples()), # bismark_mapping pe report
    output:
        file = expand("{rep_dir}/summary/03_bismark_mapping_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/03_bis_map_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/03_bis_map_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "03_bismark_mapping_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for reads mapped by bismark. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "03 Bismark Mapping QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

# 04 bwameth mapped: post deduplication
rule multiqc_compile_reports_04_bwa:
    input:
        expand("{rep_dir}/04_sambamba_dedup/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}-{sample.accession}.log", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # sambamba dedup report log for bwa mapped reads
        expand("{rep_dir}/04_fastqc_post_dedup_bwa/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bwa output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bwa/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_trimmed_sorted_dedup.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bwameth mapped
    output:
        file = expand("{rep_dir}/summary/04_bwa_deduped_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/04_bwa_deduped_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/04_bwa_deduped_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "04_bwa_deduped_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for sambamba deduplication of reads mapped by bwameth. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "04 Sambamba Deduplication (bwameth mapped) QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

# 04 bismark mapped: post deduplication
rule multiqc_compile_reports_04_bis:
    input:
        expand("{rep_dir}/04_fastqc_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated_fastqc.{suf}", rep_dir = config["reports_dir"], suf=["html","zip"], sample=sample_info.itertuples()), # fastqc bis output post-dedup
        expand("{rep_dir}/04_samtools_post_dedup_bis/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}_bismark.deduplicated.bam.stats", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # samtools_stats output bismark mapped
        expand("{rep_dir}/04_bismark_deduplication/{sample.ref}/{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}/{sample.accession}__bismark.deduplication_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]), # bismark_deduplication output report for bismark mapped reads
    output:
        file = expand("{rep_dir}/summary/04_bis_deduped_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/04_bis_deduped_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/04_bis_deduped_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "04_bis_deduped_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for Bismark deduplication of reads mapped by bismark. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "04 Bismark Deduplication (bismark mapped) QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

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
        file = expand("{rep_dir}/summary/05_bwa_post_merge_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/05_bwa_post_merge_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/05_bwa_post_merge_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "05_bwa_post_merge_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report after sambamba merging of reads mapped by bwameth. Fastqc, featureCounts, qualimap, samtools_stats, and mosdepth reports are included. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "05 Sambamba Merged (bwameth mapped) QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

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
        file = expand("{rep_dir}/summary/05_bis_post_merge_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/05_bis_post_merge_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/05_bis_post_merge_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "05_bis_post_merge_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report after sambamba merging of reads mapped by bismark. Fastqc, featureCounts, qualimap, samtools_stats, and mosdepth reports are included. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "04 Sambamba Merged (bismark mapped) QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """

# 06 bismark mapped: methylation extraction reports
rule multiqc_compile_reports_06_bis:
    input:
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.M-bias.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}_splitting_report.txt", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bismark.cov.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"]),
        expand("{rep_dir}/06_bismark_methyl_extractor/{sample.ref}-{sample.patient_id}-{sample.group}-{sample.srx_id}-{sample.layout}.bedGraph.gz", sample=sample_info.itertuples(), rep_dir=config["reports_dir"])
    output:
        file = expand("{rep_dir}/summary/06_bis_methyl_extract_multiqc.html", rep_dir=config["reports_dir"]),
        input_list = expand("{rep_dir}/summary/06_bis_methyl_extract_qc_report_list.txt", rep_dir=config["reports_dir"])
    log:
        "../pre-processing/logs/rule-logs/06_bis_methyl_extract_multiqc.log"
    conda:
        "../environment_files/multiqc.yaml"
    params:
        output_filename = "06_bis_methyl_extract_multiqc", # Required: do not change without adjusting the rule output
        outdir = config["reports_dir"] + "/summary", # Required: do not change without adjusting the rule output
        comment = "Multiqc report for bismark methylation extraction from reads mapped by bismark. Report compiled by snakemake DMR_workflow pipeline.", # Optional: comment for multiqc report (can be changed freely)
        report_title = "06 Bismark Methylation Extractor (bismark mapped) QC Reports", # Optional: title for multiqc report (can be changed freely)
        extra = "--verbose --fullnames --dirs --dirs-depth 3"  # Optional: extra parameters for multiqc (can be changed freely)
    shell:
        """
        echo "Running multiqc for raw sequence qc reports: {input}" > {log}
        echo "Output will be written to {params.outdir}/{params.output_filename}.html" >> {log}
        
        echo "create directory {params.outdir}" >> {log}
        mkdir -p {params.outdir}
        
        echo "if {output.input_list} already exists, remove it" >> {log}
        if [ -f {output.input_list} ]; then 
            rm {output.input_list}
        fi

        echo "write input list to {output.input_list}" >> {log}
        touch {output.input_list}
        for i in {input}; do 
            echo $i >> {output.input_list}
        done

        echo "run multiqc" >> {log}
        multiqc --force --filename {params.output_filename} --outdir {params.outdir} --file-list {output.input_list} --title {params.report_title} --comment {params.comment} {params.extra} >> {log} 2>&1

        echo "multiqc completed" >> {log}  
        """
