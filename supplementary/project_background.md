## Identifying differentially methylated regions (DMRs) associated with colorectal cancer

The thesis project associated with this workflow seeks to aggregate genome wide methylation data for colorectal cancer patients and healthy patients to identify common differentially methylated regions (DMRs) between sample groups.

### Background

#### DNA methylation (DNAm)
In mammals, DNAm most commonly occurs as CpG sites, where cytosine is followed immediatedly by guanine in the 5 prime to 3 prime direction.

DNAm plays a key role in determining gene expression by respressing transcriptional expression. DNAm has been shown to drive cellular identity allowing for cell differentiation. The overall levels of DNAm has been found to change with age, and abnormal DNAm patterns have been associated with various conditions including cancers.

#### Methylation patterns observed in cancers

Epigenome-wide association studies (EWAS) investigate relationships between epigenetic modifications across the entire genome and a particular condition. 

DMRs associated with a condition are identified by comparing DNA methylation (DNAm) in two groups.

<img src="https://hackmd.io/_uploads/r19wRN8UT.png" width=80% height=80%>

Global loss of DNAm is a common feature of cancer cells, with hypomethylation of intergenic repeats and gene promoters. This hypomethylation can lead to genomic instability.

Hypermethylation of gene promoters can lead to gene silencing, while hypermethylation of enhancers can reduce gene transcription. CTCF binding sites can be affected by both hypermethylation and hypomethylation, leading to genomic instability.

DNAm patterns observed in cancer cells:

| Genomic feature    | Change           | Impact of change           |
| ------------------ | ---------------- | -------------------------- |
| Intergenic repeats | Hypomethylation  | Genomic instability        |
| Gene promoters     | Hypomethylation  | Gene reactivation          |
| Gene promoters     | Hypermethylation | Gene silencing             |
| Enhancer           | Hypermethylation | Reduces gene transcription |
| CTCF binding sites | Both             | Genomic instability        |

By identifying DMRs associated with cancer, we can gain insight into the underlying mechanisms of cancer development and progression. In addition to this, we can identify targets for therapeutic intervention and biomarkers for early detection.

#### Challenges in DNAm analysis

Methylation arrays, despite covering ~2% of potential methylated regions, are commonly used due to cost considerations. 

Whole genome bisulfite sequencing (WGBS), is limited in usage due to high costs, often restricting analysis to one cancer tissue sample per study.

While comparing WGBS data between studies raises many challenges, there is value in aggregating WGBS data from various studies and analysing as a collective.

This study aims to collect open source WGBS data from colorectal cancer studies and identify common DMRs between sample groups across studies.

<!-- #### WGBS Assay
Whole genome bisulfite sequencing (WGBS) is used to investigate genome wide DNA methylation (DNAm) at a single base resolution. Isolated DNA is treated with bisulfite, converting unmethylated cytosines to uracils, while methylated cytosines remain the same. After PCR amplification and sequencing of bisulfite treated DNA, uracils will read as thymines. -->