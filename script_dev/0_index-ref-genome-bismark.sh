#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J index                                # job name
#SBATCH -t 0-5:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=15gb                              # memory pool to all cores
#SBATCH -e slurm.index.bismark.j%j.err          # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.index.bismark.j%j.out          # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

'''
------------------------
 description and usage
------------------------
0_index-ref-genome-bismark.sh
Created: 1/24/24 MS

this script will index a reference genome for alignment with bismark

input directory containing reference genome in fasta.gz format:
  - example input:
    - /home/user/genomes/human/hg38/GRCh38_no_alt_analysis_set/
    - specified location must contain reference fasta file (ex: "GRCh38_no_alt_analysis_set.fa.gz")

output is a directory of index files:
  - bismark default output directory is named "Bisulfite_Genome"
  - 3 directories, 14 files
  - structure:
        Bisulfite_Genome/
        ├── CT_conversion
        │   ├── BS_CT.1.bt2
        │   ├── BS_CT.2.bt2
        │   ├── BS_CT.3.bt2
        │   ├── BS_CT.4.bt2
        │   ├── BS_CT.rev.1.bt2
        │   ├── BS_CT.rev.2.bt2
        │   └── genome_mfa.CT_conversion.fa
        └── GA_conversion
            ├── BS_GA.1.bt2
            ├── BS_GA.2.bt2
            ├── BS_GA.3.bt2
            ├── BS_GA.4.bt2
            ├── BS_GA.rev.1.bt2
            ├── BS_GA.rev.2.bt2
            └── genome_mfa.GA_conversion.fa

options:
  -r indicates the directory containing reference genome to use (use full path to directory)
  - does not take positional arguments

usage:
sbatch 0_index-ref-genome-bismark.sh -r [reference-genome]

resources used:
indexing human hg38 analysis set (GRCh38_no_alt_analysis_set.fa.gz)
  - time: 3 hours
  - memory: 11 gb
-------------------------------------------
'''

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate bismark

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:r:' flag
do
    case "${flag}" in
        r) reference=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# Index specified file as reference genome for bismarkuse
bismark_genome_preparation --verbose --bowtie2 --parallel 6 $reference

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM

