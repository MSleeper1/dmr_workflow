# Info:
    # Script name: 0_create_annotation_references_simplified.ipynb
    # Purpose: Create annotation references for DMR analysis
    # Author: Meghan M. Sleeper

# Inputs: 
    # gene annotations
        # supplementary/gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3

    # regulatory build annotations 
        # supplementary/GRCh38.p14_regulatory_features_mart_export.tsv 

    # CGI annotations (have not added CGI annotations  to script yet)
        # supplementary/UCSC.CGI.GRCh38.bed 
    
# Outputs: 
    ### only gene annotations
        # gencode.v44.annotation.genes.tsv
        # gencode.v44.annotation.genes.sorted.formatted.bed.gz
        # gencode.v44.annotation.genes.sorted.formatted.bed.gz.tbi

    ### annotations for all features
        # gencode.v44.annotation.all.tsv (removed in last step to save space since file is 1.1 gb)
        # gencode.v44.annotation.all.sorted.formatted.bed.gz
        # gencode.v44.annotation.all.sorted.formatted.bed.gz.tbi

    ### regulatory build annotations
        # ensembl.GRch38.p14.regulatory.tsv
        # ensembl.GRch38.p14.regulatory.sorted.formatted.bed.gz
        # ensembl.GRch38.p14.regulatory.sorted.formatted.bed.gz.tbi

# Notes:
    # Modified from: https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2

# Reminders to self:
    # gencode_exons_codons_utrs_cds.head(10)
    # gencode.feature.unique()
    # gencode.head()
    # gencode[gencode.feature == 'gene'].head(10)
    # list(gencode.attribute)[0:10]
    # gencode[gencode.attribute.str.contains('gene_type')].head(10)
    # magic commands: 
        # %reset can remove variables from memory

# %%
# Import packages
import pandas as pd
import numpy as np
import pysam
import gzip
import shutil
import os

# %%
# Assign file path variables

### input file names and paths
# zipped gttf annotation reference
gencode_annotation_file = "/Users/meghansleeper/Python/deconvolution_projects/dmr_workflow/supplementary/gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3.gz"

# tab separated file with header containing regulatory build information (not sorted ot indexed yet)
ensembl_regulatory_build = "/Users/meghansleeper/Python/deconvolution_projects/dmr_workflow/supplementary/GRch38.p14_regulatory_features_mart_export.tsv"

# tab separated file with header containing CGI information (not sorted ot indexed yet)
# cgi_annotation_file = "/Users/meghansleeper/Python/deconvolution_projects/dmr_workflow/supplementary/UCSC.CGI.GRCh38.bed"

### output file names
outname_all = 'gencode.v44.annotation.all.tsv' # all features
outname_genes = 'gencode.v44.annotation.genes.tsv' # genes only
outname_reg = 'ensembl.GRch38.p14.regulatory.tsv' # regulatory build only

### output directory path
outdir = '../supplementary/annotation_beds/'

### output file names with path
fullname_all = os.path.join(outdir, outname_all) # all features
fullname_genes = os.path.join(outdir, outname_genes) # genes only
fullname_reg = os.path.join(outdir, outname_reg) # regulatory build only


# %%
#w unzip gencode_annotation_file if it is a .gz file
if gencode_annotation_file.endswith(".gz"):
    with gzip.open(gencode_annotation_file, 'rb') as f_in:
        with open(gencode_annotation_file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            gencode_annotation_file = gencode_annotation_file[:-3]

# %% [markdown]
# ## Reading in and cleaning regulatory build file

# %%
# regulatory build file to df
# note: this file only needs to be sorted and indexed, which will happen at the end of this script
annot_regulatory = pd.read_csv(ensembl_regulatory_build, 
                                  sep='\t', 
                                  header=0, names=['chr', 'start', 'end', 
                                                   'feature_type', 'feature_type_description'])


# %%
### filter out chromosomes and contigs that are not in the chromosomes list

# annot_regulatory.chr.unique() can be used to view all chromosomes in the annot_regulatory dataframe and returns the following:
        # array(['18', '8', '6', '3', '16', 'X', '7', '20', '15', '19', '5', '2',
        #        '22', '21', '17', '1', '4', '10', '12', '14', '9', '11', '13', 'Y',
        #        'GL000195.1', 'KI270731.1', 'KI270736.1', 'KI270712.1',
        #        'GL000214.1', 'GL000194.1', 'KI270729.1', 'KI270333.1',
        #        'KI270728.1', 'GL000224.1', 'KI270742.1', 'KI270744.1',
        #        'GL000220.1', 'KI270737.1', 'KI270727.1', 'KI270706.1',
        #        'KI270442.1', 'GL000221.1', 'GL000219.1', 'KI270721.1',
        #        'KI270725.1', 'KI270438.1', 'KI270752.1', 'KI270756.1',
        #        'KI270730.1', 'KI270734.1', 'KI270716.1', 'KI270711.1',
        #        'GL000216.2', 'KI270719.1', 'KI270739.1', 'KI270751.1',
        #        'GL000205.2', 'KI270743.1', 'KI270747.1', 'KI270738.1',
        #        'GL000218.1', 'GL000208.1', 'KI270733.1', 'KI270717.1',
        #        'KI270741.1', 'GL000008.2', 'GL000213.1', 'KI270757.1',
        #        'KI270722.1', 'GL000009.2', 'KI270750.1', 'KI270748.1',
        #        'KI270749.1', 'KI270745.1', 'KI270519.1', 'KI270710.1',
        #        'KI270754.1', 'GL000225.1', 'KI270713.1', 'KI270714.1',
        #        'KI270528.1', 'KI270707.1', 'KI270723.1', 'KI270718.1',
        #        'KI270732.1', 'KI270720.1', 'KI270753.1', 'KI270708.1',
        #        'KI270521.1', 'KI270709.1', 'KI270538.1', 'KI270735.1',
        #        'KI270581.1', 'KI270726.1', 'KI270755.1', 'KI270724.1',
        #        'KI270336.1', 'KI270330.1', 'KI270466.1', 'KI270337.1'],
        #       dtype=object)

# list of chromosome numbers to include in final annotation file (excludes contigs and alternative haplotypes)
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 
               '10', '11', '12', '13', '14', '15', '16', '17', 
               '18', '19', '20', '21', '22', 'X', 'Y']

# if annot_regulatory['chr'] is not in chromosomes list, remove the row
annot_regulatory = annot_regulatory[annot_regulatory['chr'].isin(chromosomes)]

# add "chr" string to the value to match the format of wgbstools output
annot_regulatory['chr'] = 'chr' + annot_regulatory['chr']


# %% [markdown]
# ## Preparing gencode gtf file for use in annotation script

# %%
# Create dataframe from gencode gtf annotation file
gencode = pd.read_table(gencode_annotation_file, 
                        comment="#", 
                        sep = "\t", 
                        names = ['seqname', 'source', 'feature', 
                                'start' , 'end', 'score', 'strand', 
                                'frame', 'attribute'
                                ]
                        )

# %%
# creating dataframes for grouped features that share the same attribute information to make parsing attributes easier
  # see https://www.gencodegenes.org/pages/data_format.html for attributes associated with each feature
    # gencode.feature.unique() returns:
        #    'gene', 'transcript', 'exon', 'CDS', 
        #    'start_codon', 'stop_codon',
        #    'five_prime_UTR', 'three_prime_UTR',
        #    'stop_codon_redefined_as_selenocysteine'

gencode_genes = gencode[(gencode.feature == "gene")][
    ['seqname', 'start', 'end', 'attribute','feature', 
     'strand', 'source']].copy().reset_index().drop('index', axis=1)

gencode_exons_codons_utrs_cds = gencode[
    (gencode.feature == "exon") | 
    (gencode.feature == "start_codon") | 
    (gencode.feature == "stop_codon") | 
    (gencode.feature == "five_prime_UTR") | 
    (gencode.feature == "three_prime_UTR") | 
    (gencode.feature == "CDS")
    ][['seqname', 'start', 'end', 'attribute', 'feature', 'strand', 'source']
      ].copy().reset_index().drop('index', axis=1)

gencode_transcripts_selenocyteines = gencode[
    (gencode.feature == "transcript") | 
    (gencode.feature == "stop_codon_redefined_as_selenocysteine")
    ][['seqname', 'start', 'end', 'attribute', 'feature', 'strand', 'source']
      ].copy().reset_index().drop('index', axis=1)


# %%
### functions for parsing gtf files based on format of attribute column
# list(gencode.attribute)[100:150] can be used to view the attribute format in the gencode dataframe

def gene_info(x):
# Extract gene names, gene id gene_type, level, and gene_id for each gene
# Use this function when targetting info on genes
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
    g_level = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
    return (g_name, g_id, g_type, g_level)

def all_info(x):
# Extract gene names, gene id, gene_type, level, exon_id, exon_number, transcript_name, and transcript_type for each exon
# Use this function when targetting info on exons, CDS, start_codon, stop_codon, five_prime_UTR, and three_prime_UTR
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
    g_level = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
    exon_id = list(filter(lambda x: 'exon_id' in x,  x.split(";")))[0].split("=")[1]
    exon_number = int(list(filter(lambda x: 'exon_number' in x,  x.split(";")))[0].split("=")[1])
    transcript_name = list(filter(lambda x: 'transcript_name' in x,  x.split(";")))[0].split("=")[1]
    transcript_type = list(filter(lambda x: 'transcript_type' in x,  x.split(";")))[0].split("=")[1]
    return (g_name, g_id, g_type, g_level, exon_id, exon_number, transcript_name, transcript_type)
    
def transcript_info(x):
# Extract gene names, gene id, gene_type, level, transcript_id, and transcript_type for each transcript
# Use this function when targetting info on transcripts and stop_codon_redefined_as_selenocysteine
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
    g_level = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
    transcript_name = list(filter(lambda x: 'transcript_name' in x,  x.split(";")))[0].split("=")[1]
    transcript_type = list(filter(lambda x: 'transcript_type' in x,  x.split(";")))[0].split("=")[1]
    return (g_name, g_id, g_type, g_level, transcript_name, transcript_type)

# %%
## Using gene_info function on genes
# extract gene_name, gene_id, gene_type, and gene_level from gencode_genes.attribute
gencode_genes["gene_name"], gencode_genes["gene_id"], gencode_genes["gene_type"], gencode_genes["gene_level"] = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))

## Using all_info function on exons, start codons, stop codons, utrs, and cds
# extract gene_name, gene_type, gene_level, exon_id, exon_number, transcript_name, and transcript_type from gencode_exons_codons_utrs_cds.attribute
gencode_exons_codons_utrs_cds["gene_name"], gencode_exons_codons_utrs_cds["gene_id"], gencode_exons_codons_utrs_cds["gene_type"], gencode_exons_codons_utrs_cds["gene_level"], gencode_exons_codons_utrs_cds["exon_id"], gencode_exons_codons_utrs_cds["exon_number"], gencode_exons_codons_utrs_cds["transcript_name"], gencode_exons_codons_utrs_cds["transcript_type"] = zip(*gencode_exons_codons_utrs_cds.attribute.apply(lambda x: all_info(x)))

## Using transcript_info function on transcripts and stop_codon_redefined_as_selenocysteine
# extract gene_name, gene_type, gene_level, exon_id, exon_number, transcript_name, and transcript_type from gencode_transcripts_selenocyteines.attribute using transcript_info function
gencode_transcripts_selenocyteines["gene_name"], gencode_transcripts_selenocyteines["gene_id"], gencode_transcripts_selenocyteines["gene_type"], gencode_transcripts_selenocyteines["gene_level"], gencode_transcripts_selenocyteines["transcript_name"], gencode_transcripts_selenocyteines["transcript_type"] = zip(*gencode_transcripts_selenocyteines.attribute.apply(lambda x: transcript_info(x)))

# %%
### genes can also be filtered by gene_type if interested in only protein coding genes

# gencode_genes.gene_type.unique() returns all of the gene types in the gencode_genes dataframe and returns the following:
    # array(['lncRNA', 'transcribed_unprocessed_pseudogene',
    #        'unprocessed_pseudogene', 'miRNA', 'protein_coding',
    #        'processed_pseudogene', 'snRNA',
    #        'transcribed_processed_pseudogene', 'misc_RNA', 'TEC',
    #        'transcribed_unitary_pseudogene', 'snoRNA', 'scaRNA',
    #        'rRNA_pseudogene', 'unitary_pseudogene', 'pseudogene', 'rRNA',
    #        'IG_V_pseudogene', 'scRNA', 'IG_V_gene', 'IG_C_gene', 'IG_J_gene',
    #        'sRNA', 'ribozyme', 'translated_processed_pseudogene', 'vault_RNA',
    #        'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_V_pseudogene',
    #        'TR_D_gene', 'IG_C_pseudogene', 'TR_J_pseudogene',
    #        'IG_J_pseudogene', 'IG_D_gene', 'IG_pseudogene', 'artifact',
    #        'Mt_tRNA', 'Mt_rRNA'], dtype=object)

# For filtering out the genes that are protein coding into a new dataframe, uncomment the following line
# gencode_genes_protein_coding = gencode_genes[gencode_genes['gene_type'] == 'protein_coding'].reset_index().drop('index', axis=1)

# %%

# concatenating dataframes for all features into one dataframe with NaNs for missing values where features don't have the same attributes
gencode_all = pd.concat(
    [gencode_genes, 
     gencode_exons_codons_utrs_cds, 
     gencode_transcripts_selenocyteines], 
     ignore_index=True
     ).reset_index().drop('index', axis=1)

# To confirm that all features are included in the new dataframe, uncomment the following line
# gencode_all.feature.unique()


# %% [markdown]
# ## Saving, sorting, and indexing new annotation reference files
# 
# The gene annotation files will be saved in bed format with tab delimition and no header row. The rows will be sorted by chromosome, start, and end (in that order).
# 
# `gencode_all` contains all feature types and will be saved  as `gencode.v44.annotation.all.sorted.formatted.bed` with the columns:
#     ```['seqname', 'start', 'end', 'attribute', 'feature', 'strand', 'source', 'gene_name','gene_id', 'gene_type', 'gene_level', 'exon_id', 'exon_number', 'transcript_name', 'transcript_type']```
# 
# `gencode_gene` contains genes only and will be saved  as `gencode.v44.annotation.genes.sorted.formatted.bed` with the columns:
#     ```['seqname', 'start', 'end', 'attribute', 'feature', 'strand', 'source', 'gene_name','gene_id', 'gene_type', 'gene_level']```
# 
# `annot_regulatory` contains regulatory features/associated loci and will be saved as `ensembl.GRch38.p14.regulatory.sorted.formatted.bed` with the columns:
#     ```['chr', 'start', 'end', 'feature_type', 'feature_type_description']```

# to find the column names for a dataframe, uncomment the following line and change dataframe name
# list(gencode_all.columns)


# %%
### Save the gencode dataframes into files to be sorted and converted to bed files

# create output directories if they don't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

#### Save the gencode dataframes into files to be sorted and converted to bed files
gencode_all.to_csv(fullname_all, index=False, header = False, sep="\t")
gencode_genes.to_csv(fullname_genes, index=False, header = False, sep="\t")
annot_regulatory.to_csv(fullname_reg, index=False, header = False, sep="\t")

# %% [markdown]

# ## Sorting, compressing, and indexing annotation reference files

# %% 
%%bash -s ../supplementary/annotation_beds/ensembl.GRch38.p14.regulatory.tsv
# Sort by chromosome, start, and end
cut -f 1,2,3,4,5 $1 | sort -k1,1 -k2,2n -k3,3n > ../supplementary/annotation_beds/ensembl.GRch38.p14.regulatory.sorted.formatted.bed

# Compress the bed file
bgzip -f ../supplementary/annotation_beds/ensembl.GRch38.p14.regulatory.sorted.formatted.bed

# Index the bed file
tabix -f -p bed ../supplementary/annotation_beds/ensembl.GRch38.p14.regulatory.sorted.formatted.bed.gz

# %%
%%bash -s ../supplementary/annotation_beds/gencode.v44.annotation.all.tsv
# Sort by chromosome, start, and end
cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 $1 | sort -k1,1 -k2,2n -k3,3n > ../supplementary/annotation_beds/gencode.v44.annotation.all.sorted.formatted.bed

# Compress the bed file
bgzip -f ../supplementary/annotation_beds/gencode.v44.annotation.all.sorted.formatted.bed

# Index the bed file
tabix -f -p bed ../supplementary/annotation_beds/gencode.v44.annotation.all.sorted.formatted.bed.gz

# %%
%%bash -s ../supplementary/annotation_beds/gencode.v44.annotation.genes.tsv
# Sort by chromosome, start, and end
cut -f 1,2,3,4,5,6,7,8,9,10,11 $1 | sort -k1,1 -k2,2n -k3,3n > ../supplementary/annotation_beds/gencode.v44.annotation.genes.sorted.formatted.bed

# Compress the bed file
bgzip -f ../supplementary/annotation_beds/gencode.v44.annotation.genes.sorted.formatted.bed

# Index the bed file
tabix -f -p bed ../supplementary/annotation_beds/gencode.v44.annotation.genes.sorted.formatted.bed.gz

# %%
### remove large intermediate files to save space

# remove unzipped gencode annotation file to save space
if not gencode_annotation_file.endswith(".gz"):
    os.remove(gencode_annotation_file)

# tsv with all features is 1.1 gb and should be removed to save space
# remove large outname_all file to save space if it was created
if os.path.exists(fullname_all):
    os.remove(fullname_all)

# %%
